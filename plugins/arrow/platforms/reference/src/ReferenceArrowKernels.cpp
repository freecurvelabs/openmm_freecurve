/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "ReferenceArrowKernels.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceConstraints.h"
#include "ReferenceVirtualSites.h"
#include <set>

 // Arbalest includes:

#include "ParseCommandLine.h"
#include "SystemLoadingLib/SystemLoader.h"
#include "SimulationCoreLib/SimulationEnvironment.h"
#include "SimulationCoreLib/SimulationEnvironmentGpu.h"
#include "SimulationCoreLib/MDSchemeOperations.h"
#include "SimulationToolsLib/SimulationController.h"
#include "GpuLib/GpuManager.h"
#include "CommonLib/Timer.h"
#include "CommonLib/Typedefs.h"
#include "CommonLib/SubprocessRun.h"
#include "SimulationCoreLib/CalculatedValues.h"

using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->positions;
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->velocities;
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->forces;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->constraints;
}

static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& inverseMasses, double timeShift) {
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    vector<Vec3> shiftedVel(numParticles);
    context.computeShiftedVelocities(timeShift, shiftedVel);
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (inverseMasses[i] > 0)
            energy += (shiftedVel[i].dot(shiftedVel[i]))/inverseMasses[i];
    return 0.5*energy;
}

bool ReferenceCalcArrowForceKernel::copyCrdFromContextToArbalest(ContextImpl& context, SimulationCore::CEnvironmentReplica* pEnvReplica)
{
    std::vector<OpenMM::Vec3> positions;
    OpenMM::Vec3 a, b, c; //Periodic box vectors

    context.getPositions(positions);
    context.getPeriodicBoxVectors(a, b, c);
    int natoms = positions.size();

    for(int i = 0; i < natoms; i++)
      positions[i] = positions[i]* 10.0; 
    a = a * 10.0;
    b = b * 10.0;
    c = c * 10.0;

    bool bRes = true;
    try
    {
        if (bRes)
        {
            double fLX, fLY, fLZ;
            pEnvReplica->m_pDomainController->GetDomainSize(fLX, fLY, fLZ);
            //We should not expect the coincidence: Barostat etc. can slightly change the domain dimensions. 
            if (fabs(a[0] - fLX) > ((double)0.025) * (fabs(a[0]) + fabs(fLX)) //Allow relative difference up to 2.5+2.5 = 5%
                || fabs(b[1] - fLY) > ((double)0.025) * (fabs(b[1]) + fabs(fLY))
                || fabs(c[2] - fLZ) > ((double)0.025) * (fabs(c[2]) + fabs(fLZ))
                )
            { // Show warning message
                TRACE_DEBUG2(_ComStr("Atom Positions Loading"), BLDSTRING13("Box dimensions in OpenMM context are different from that in the system: (", a[0], ", ", b[1], ", ", c[2], ") vs. (", fLX, ", ", fLY, ", ", fLZ, ") in Config"));
            }
            //[Slightly] correct box dimensions with values loaded from TRR
            {
                double fScaleFactorX = a[0] / fLX;
                double fScaleFactorY = b[1] / fLY;
                double fScaleFactorZ = c[2] / fLZ;
                bRes = pEnvReplica->m_pDomainController->RescaleDomainSize(fScaleFactorX, fScaleFactorY, fScaleFactorZ); //Be careful: this does not recalculate hash table
                pEnvReplica->m_pDomainController->GetDomainSize(fLX, fLY, fLZ); // Запрашиваем [обновленные] размеры бокса
                TRACE_DEBUG2(_ComStr("Atom Positions Loading"), BLDSTRING7("Box dimensions were replaced by values of (", fLX, ", ", fLY, ", ", fLZ, ") taken from OpenMM context"));
            }
        }
        if (bRes && natoms != pEnvReplica->m_pCmptAtoms->m_nAtoms)
        {
            TRACE_FATAL2(_ComStr("Atom Positions Loading"), BLDSTRING5("Number of atoms in OpenMM context (", natoms, ") is different from that in the system (", pEnvReplica->m_pCmptAtoms->m_nAtoms, ")"));
            bRes = false;
        }

        if (bRes)
        {
            // False => Atoms are supposed to be stored in the 'original' order, as they initially came to the system from the structure molecular files
            // True => Atoms are supposed to be stored in the sorted order. Hopefully that sorting order was the same as current order in CmptAtoms, i.e. 'internal' order. So just copy arrays 'as is'
            bool bSortedAtomsOrder = true;
            SimulationCore::CFunctionalGroup* pFuncGroupSystem = pEnvReplica->m_pFuncGroups->GetGroupSystem();
            Common::IndexType iAtom = -1, nUsedAtomsNumber = 0;
            const Common::IndexType* pUsedAtomicIndices = pFuncGroupSystem->GetAtomicIndices(bSortedAtomsOrder, nUsedAtomsNumber);
            for (int i = 0; i < natoms; i++)
            {
                //if( i < 3 ) {
                //   std::cout << " Atoms Positions: " << std::endl; 
                //   std::cout << i << "  " << positions[i][0]  << "  " << positions[i][1] << "  " << positions[i][2] << std::endl;
                //}
                pEnvReplica->m_pCmptAtoms->m_pR[i].x = DBL_TO_COORD( positions[i][0] );
                pEnvReplica->m_pCmptAtoms->m_pR[i].y = DBL_TO_COORD( positions[i][1] );
                pEnvReplica->m_pCmptAtoms->m_pR[i].z = DBL_TO_COORD( positions[i][2] );
            }
        }
    }
    catch (...)
    {
        TRACE_FATAL2(_ComStr("Atom Positions Loading"), BLDSTRING1("Unknown exception while loading atoms positions from OPENMM context "));
        bRes = false;
    }

    return bRes;
}


void ReferenceCalcArrowForceKernel::initialize(const System& system, const ArrowForce& force) {
    // Initialize particle parameters.
    
    std::cout << "ReferenceCalcArrowForceKernel::initialize() pt 1 " << std::endl;
    this->scale_force = force.scale_force;
    std::cout << "Initialize Arbalest structures with config file: " << force.getArbalestConfig() << std::endl;

    int argc = 5;
    char** argv = new char* [argc];

    char argv_2[256];

    argv[0] = "Arbalest_openmm";
    argv[1] = "--config";
    strcpy(argv_2, force.getArbalestConfig().c_str());
    argv[2] = argv_2;
    argv[3] = "--gpu";  // For GPU support - to move to CudaArrowKernel
    argv[4] = "1";

    int nMyID = -1, nNumProcs = -1;
    std::string strProcTitle;
    // MPI
#ifdef ARBALEST_MPI
    int  nTitleLen;
    char sProcTitle[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nNumProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &nMyID);
    MPI_Get_processor_name(sProcTitle, &nTitleLen);
    strProcTitle = sProcTitle;
#endif	// MPI		

    try
    {
        Arbalest::CmdLineContext CmdParams;
        //        if (Arbalest::parse_command_line(argc, argv, CmdParams))
        //        {
        //            std::string sCurConfigFile = CmdParams.sConfigFile;
        //            cout << "Arbalest Config file " << sCurConfigFile << std::endl;
        //        } 
        bool useGpu = false;
        //bool useGpu = true; 

        pSysLdr = std::make_shared< SystemLoading::CSystemLoader >(useGpu);
        // SystemLoading::CSystemLoader SysLdr(useGpu);
        std::string sCompilationDetails = "";
        CmdParams.iOpenMPThreads = 8;
        CmdParams.bDumpConf = true;
        CmdParams.bMkDirs = true;
        CmdParams.bOutputNNDescriptorsToFile = false;
        CmdParams.nOutputNNBesselDescriptorsToFile = false;

        if (!pSysLdr->LoadSystem(force.getArbalestConfig(), sCompilationDetails, CmdParams.iOpenMPThreads, CmdParams.bDumpConf, CmdParams.bMkDirs
            , CmdParams.bOutputNNDescriptorsToFile // For NN
            , CmdParams.nOutputNNBesselDescriptorsToFile // For NNBessel
        ))
        {
            // LOG_INFO << "Arbalest failed to load system configuration!" << FLUSH_LOG;
            std::cout << "Arbalest failed to load system configuration!" << std::endl;
            //gpuinfoOut(GPUTestInfoOutput::eTest).OutLn("    <failure>Failed to load config file</failure>");
        }
        else
        {
            // LOG_INFO << "Simulation configuration successfully loaded!" << FLUSH_LOG;
            std::cout << "Simulation configuration successfully loaded!" << std::endl;
        }

        if (!SimController.Initialize(pSysLdr->GetSimulationReferences(), pSysLdr->GetSimulationenvironment(), pSysLdr->GetTaskContainer(), CmdParams.bOutTimeStamp))
        {
            LOG_INFO << "Simulation initialization failed!" << FLUSH_LOG;
            // gpuinfoOut(GPUTestInfoOutput::eTest).OutLn("    <failure>Simulation init failed</failure>");
            // return 0;
            return;
        }
        if (CmdParams.bVerify)
        {
#ifdef ARBALEST_CUDA
            pSysLdr->GetSimulationenvironment()->EnableVerificationGPUvsCPU();
#endif // ARBALEST_CUDA
        }
        if (CmdParams.bGPUSynchronize)
        {
            SimulationCore::CSimulationEnvironmentGPU* pSimEnvGpu = dynamic_cast<SimulationCore::CSimulationEnvironmentGPU*>(pSysLdr->GetSimulationenvironment());
            if (pSimEnvGpu)
            {
                pSimEnvGpu->EnableSynchronizationGPUvsCPU(true);
            }
        }

        SimController.m_pTaskContainer->Initialize(SimController.m_pSimRefs, SimController.m_pSimEnv);
		SimController.PreLaunchSimulation(); // Moving Simulation setup here    
		
		// Move functions from CEnergyValuation::Launch()
		// Notify objects about beginning of the task

        std::cout << "ReferenceCalcArrowForceKernel::initialize() before SimController.m_pSimEnv->OnTaskStarted() pt 10 " << std::endl; 
		SimController.m_pSimEnv->OnTaskStarted();  
		spMDSchemeOperations = SimulationCore::CreateMDSchemeOperationsObject(SimController.m_pSimEnv);  
    }
    catch (const Common::CException& e)
    {
        // LOG_FATAL << "Common::CException: " << e.Message() << FLUSH_LOG;
        std::cout << "Common::CException: " << e.Message() << std::endl;
        //throw e; // Do not throw and return 'correctly'
    }
    catch (std::runtime_error& e)
    {
        // LOG_FATAL << "std::runtime_error: " << e.what() << FLUSH_LOG;
        std::cout << "std::runtime_error: " << e.what() << std::endl;
        //throw e; // Do not throw and return 'correctly'
    }

    //cc.addForce(new CommonArrowForceInfo(force));
    delete[] argv;
    
}

bool bFirstTimePairs = true;

double ReferenceCalcArrowForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //cout << "ReferenceCalcArrowForceKernel::execute()" << std::endl;
    SimulationCore::CSimulationEnvironment* pEnv = pSysLdr->GetSimulationenvironment();

    bool bRes = copyCrdFromContextToArbalest(context, pEnv);

    try
    {
        // bool bSimulationResult = SimController.LaunchSimulation();
		
		// Expansion of  SimController.LaunchSimulation()
		bool bUseCmptThreadAlgorithms = false;
        std::vector<std::shared_ptr<SimulationCore::CAlgorithm> >* pvCmptThAlgorithms = NULL;
        //bool bPerformTaskLaunches = true;
        //bRes = bRes && SimController.m_pTaskContainer->LaunchSequence();
		// Expand LaunchSequence for SimulationTools::CEnergyValuation::Launch()
	
         SimulationCore::ECreateAtomPairsHints eCreateAtomPairsHint = SimulationCore::eMandatoryCreateAtomPairs;  // Recompute pairs on every step?  So far this looks more stable
    //   SimulationCore::ECreateAtomPairsHints eCreateAtomPairsHint = SimulationCore::eCreateAtomPairsOnlyIfNecessary;
        bool bCalculateAggregatedValues = true;
    //    bool bCalculateAggregatedValues = false;
        SimulationCore::CBARDynamics* pBARDynamics = NULL;
        
        bRes = spMDSchemeOperations->ComputeEnergyAndForces(SimController.m_pSimEnv, bFirstTimePairs, eCreateAtomPairsHint, bCalculateAggregatedValues, pBARDynamics); 
		//bFirstTimePairs = false;   
		bFirstTimePairs = true;  // Recompute pairs on every step?  So far this looks more stable
    }
    catch (const Common::CException& e)
    {
        // LOG_FATAL << "Common::CException: " << e.Message() << FLUSH_LOG;
        std::cout << "Common::CException: " << e.Message() << std::endl;
        //throw e; // Do not throw and return 'correctly'
    }
    catch (std::runtime_error& e)
    {
        // LOG_FATAL << "std::runtime_error: " << e.what() << FLUSH_LOG;
        std::cout << "std::runtime_error: " << e.what() << std::endl;
        //throw e; // Do not throw and return 'correctly'
    }

    SimulationCore::FORCEVEC3* pFrc = pEnv->m_pCmptAtoms->m_pFtotal;

    //printf(" ReferenceCalcArrowForceKernel::execute() line 303  pFrc[0].x = %f ", FORCE_TO_DBL(pFrc[0].x));
    //printf("pFrc[0].y = %f ", FORCE_TO_DBL(pFrc[0].y)); 
    //printf("pFrc[0].z = %f \n", FORCE_TO_DBL(pFrc[0].z));    

    double energy = 0;
    vector<Vec3>& forces = extractForces(context);
    int natoms = forces.size();

    //for( int i = 0; i < 2; i++)
    //{
    //    printf(" ReferenceCalcArrowForceKernel::execute() line 311  forces[%d] = %f %f %f \n", i, forces[i][0], forces[i][1], forces[i][2]);
    //}
    //printf(" ReferenceCalcArrowForceKernel::execute() line 315 \n" );

    for (int i = 0; i < natoms; i++)
    {
        forces[i][0] += FORCE_TO_DBL( pFrc[i].x ) * 41.84 * scale_force;
        forces[i][1] += FORCE_TO_DBL( pFrc[i].y ) * 41.84 * scale_force;
        forces[i][2] += FORCE_TO_DBL( pFrc[i].z ) * 41.84 * scale_force;
    }

    // Potential energy in kJ/mol - seems to be defaul in OpenMM ??
    energy = pEnv->m_pCmptMatrix->m_pfCalculatedValues[0][SimulationCore::ECalculatedValues::eEnergPot]*4.184*scale_force;

    // std::cout << " ReferenceCalcArrowForceKernel::execute() : Potential Energy =" <<  energy/4.184 << "  kcal/mol " << std::endl; 
	// std::cout << " ReferenceCalcArrowForceKernel::execute() : Potential Energy =" <<  energy << "  kJ/mol " << std::endl;   
    return energy;
}

void ReferenceCalcArrowForceKernel::copyParametersToContext(ContextImpl& context, const ArrowForce& force) {
    if (force.getNumParticles() != particle.size())
        throw OpenMMException("updateParametersInContext: The number of Drude particles has changed");

    for (int i = 0; i < force.getNumParticles(); i++) {
        int p;
        force.getParticleParameters(i, p);
        if (p != particle[i] )
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}

/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "ReferencePlatform.h"
#include "CudaPlatform.h"
#include "CudaContext.h"

#include "CommonArrowKernels.h"
#include "CommonArrowKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/CommonKernels.h"
#include "openmm/common/IntegrationUtilities.h"
//#include "openmm/cuda/CudaPlatform.h"
//#include "openmm/cuda/CudaContext.h"
//#include "openmm/cuda/CudaContext.h"
#include "CommonKernelSources.h"
#include "SimTKOpenMMRealType.h"
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

class CommonArrowForceInfo : public ComputeForceInfo {
public:
    CommonArrowForceInfo(const ArrowForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumParticles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        particles.clear();
        if (index < force.getNumParticles()) {
            int p;
            force.getParticleParameters(index, p);
            particles.push_back(p);
        }
    }
    bool areGroupsIdentical(int group1, int group2) {
        if (group1 < force.getNumParticles() && group2 < force.getNumParticles()) {
            return (true);
        }
        if (group1 >= force.getNumParticles() && group2 >= force.getNumParticles()) {
            return true;
        }
        return false;
    }
private:
    const ArrowForce& force;
};

bool CommonCalcArrowForceKernel::copyCrdFromContextToArbalest(ContextImpl& context, SimulationCore::CEnvironmentReplica* pEnvReplica)
{
    std::vector<OpenMM::Vec3> positions;
    OpenMM::Vec3 a, b, c; //Periodic box vectors

    context.getPositions(positions);
    context.getPeriodicBoxVectors(a, b, c);
    int natoms = positions.size();

    for(int i = 0; i < natoms; i++)
      positions[i] = positions[i] * 10.0; 

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

            SimulationCore::VECVAL3D *pvecAtomShift = &pSysLdr->GetSimulationenvironment()->m_vecAtomShift;

            for (int i = 0; i < natoms; i++)
            {
                //if( i < 3 ) {
                //   std::cout << " Atoms Positions: " << std::endl; 
                //   std::cout << i << "  " << positions[i][0]  << "  " << positions[i][1] << "  " << positions[i][2] << std::endl;
                //}
                pEnvReplica->m_pCmptAtoms->m_pR[i].x = DBL_TO_COORD( positions[i][0] - VEC_X(pvecAtomShift[0]) );
                pEnvReplica->m_pCmptAtoms->m_pR[i].y = DBL_TO_COORD( positions[i][1] - VEC_Y(pvecAtomShift[0]) );
                pEnvReplica->m_pCmptAtoms->m_pR[i].z = DBL_TO_COORD( positions[i][2] - VEC_Z(pvecAtomShift[0]) );
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

OpenMM::Vec3 CommonCalcArrowForceKernel::getArbalestShiftVec() const
{
    OpenMM::Vec3 shift; // Arbalest environment shift vector for atoms in 
    SimulationCore::VECVAL3D *pvecAtomShift = &pSysLdr->GetSimulationenvironment()->m_vecAtomShift;
    shift[0] = VEC_X(pvecAtomShift[0]) * 0.1;
    shift[1] = VEC_Y(pvecAtomShift[0]) * 0.1;
    shift[2] = VEC_Z(pvecAtomShift[0]) * 0.1; // Convert to nanometers

    return shift;
}

bool CommonCalcArrowForceKernel::saveInternalPositions( SimulationCore::CEnvironmentReplica* pEnvReplica, ContextImpl& context )
{
    // printf("CommonCalcArrowForceKernel::copyCrdFromArbalestToContext() \n");
    OpenMM::Vec3 a, b, c; //Periodic box
    std::vector<OpenMM::Vec3> positions_old;

    //context.getPeriodicBoxVectors(a, b, c);
    context.getPositions(positions_old);
    
    int natoms = pEnvReplica->m_pCmptAtoms->m_nAtoms;
    positions_internal.resize(natoms, OpenMM::Vec3(0.0, 0.0, 0.0));
    
    bool bRes = true;
    SimulationCore::VECVAL3D *pvecAtomShift = &pSysLdr->GetSimulationenvironment()->m_vecAtomShift;
    for (int i = 0; i < natoms; i++)
    {
        positions_internal[i][0] = COORD_TO_DBL(pEnvReplica->m_pCmptAtoms->m_pR[i].x) + VEC_X(pvecAtomShift[0]);
        positions_internal[i][1] = COORD_TO_DBL(pEnvReplica->m_pCmptAtoms->m_pR[i].y) + VEC_Y(pvecAtomShift[0]);
        positions_internal[i][2] = COORD_TO_DBL(pEnvReplica->m_pCmptAtoms->m_pR[i].z) + VEC_Z(pvecAtomShift[0]);
    }

    this->positions_changed = false;
    for(int i = 0; i < natoms; i++)
    {
        positions_internal[i] = positions_internal[i] * 0.1; // Convert to nanometers

        //double d2 = 0.0;
        //for (int j = 0; j < 3; j++) {
        //    d2 += (positions_internal[i][j] - positions_old[i][j]) * (positions_internal[i][j] - positions_old[i][j]);
        //}
        
        //if( d2 > 0.0001 && i < 1000) { // 0.0001 nm^2 = 0.01 Angstroms
        //    positions_changed = true;
        //    std::cout << " CommonCalcArrowForceKernel::copyCrdFromArbalestToContext()  Atoms Positions: " << std::endl; 
        //    std::cout << i << " old: " << positions_old[i][0]  << "  " << positions_old[i][1] << "  " << positions_old[i][2] << std::endl;
        //    std::cout << i << " new: " << positions_internal[i][0]  << "  " << positions_internal[i][1] << "  " << positions_internal[i][2] << std::endl;
        //}
    }

    //if( positions_changed ) {
    //    this->positions_changed = true;
    //    context.setPositions(positions);
    //}
    // context.setPeriodicBoxVectors(a, b, c); No changing of periodic box vectors so far

    return bRes;
}


void CommonCalcArrowForceKernel::initialize(const System& system, const ArrowForce& force) {
    // Initialize particle parameters.
    
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
        //bool useGpu = false;   // Do not use GPU for now
        
        std::string platform_name = kernel_platform.getName();

        bool useGpu = false; 
        if (platform_name == "CUDA") {  
            useGpu = true;
        }

        pSysLdr = std::make_shared< SystemLoading::CSystemLoader >(useGpu);  

        // SystemLoading::CSystemLoader SysLdr(useGpu);
        std::string sCompilationDetails = ""; 
        CmdParams.iOpenMPThreads = 4;  
        CmdParams.bDumpConf = true;
        CmdParams.bVerify = false;
        CmdParams.bMkDirs = true;
        CmdParams.bOutputNNDescriptorsToFile = false;
        CmdParams.nOutputNNBesselDescriptorsToFile = -1; // -1 means no output for NNBessel descriptors
        //printf("CmdParams.iGpuDevId = %d\n", CmdParams.iGpuDevId);
        //printf("CmdParams.bUseGpu = %d\n", CmdParams.bUseGpu);
        //printf("CmdParams.bOutTimeStamp = %d\n", CmdParams.bOutTimeStamp);
        //printf("CmdParams.bGPUTestInfo = %d\n", CmdParams.bGPUTestInfo);
        //printf("CmdParams.bGPUSynchronize = %d\n", CmdParams.bGPUSynchronize);
        //printf("CmdParams.iLogLevel= %d\n", CmdParams.iLogLevel);
        //printf("CmdParams.iMpiProcDelay= %d\n", CmdParams.iMpiProcDelay);
        //printf("CmdParams.sConfigFile= %s\n", CmdParams.sConfigFile.c_str());

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

        SimulationCore::VECVAL3D *pvecAtomShift = &pSysLdr->GetSimulationenvironment()->m_vecAtomShift;

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


void CommonCalcArrowForceKernel::addForces(vector<Vec3>& forces_loc, ContextImpl& context) {

    std::string platform_name = context.getPlatform().getName();
    int natoms = context.getSystem().getNumParticles();

    if (platform_name == "Reference" || platform_name == "CPU")
    {
        ReferencePlatform::PlatformData *data = reinterpret_cast<ReferencePlatform::PlatformData *>(context.getPlatformData());
        for(int i = 0; i < forces_loc.size(); i++)
        {
            (*(data->forces))[i][0] += forces_loc[i][0];
            (*(data->forces))[i][1] += forces_loc[i][1];
            (*(data->forces))[i][2] += forces_loc[i][2];
        }
    }
    else if (platform_name == "CUDA")
    {
        ContextSelector selector(*p_cc);
        long long *force = (long long *)p_cc->getPinnedBuffer();
        p_cc->getLongForceBuffer().download(force);
        const vector<int> &order = p_cc->getAtomIndex();
        int numParticles = context.getSystem().getNumParticles();
        int paddedNumParticles = p_cc->getPaddedNumAtoms();
        // forces.resize(numParticles);
        double scale = (double)0x100000000LL;
        for (int i = 0; i < numParticles; ++i)
        {
            force[i]                          += forces_loc[order[i]][0] * scale;
            force[i + paddedNumParticles]     += forces_loc[order[i]][1] * scale;
            force[i + paddedNumParticles * 2] += forces_loc[order[i]][2] * scale;
        }
        p_cc->getLongForceBuffer().upload(force);
    }
    else
    {
        printf(" CommonCalcArrowForceKernel::setForces() - Unsupported platform: %s \n", platform_name.c_str());
        throw OpenMMException("CommonCalcArrowForceKernel::execute() - Unsupported platform");
    }
}

bool bFirstTimePairs = true;

double CommonCalcArrowForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
        bRes = saveInternalPositions(pEnv, context); // save coordinates from Arbalest environment in positions_arbalest
		 
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

    // printf(" CommonCalcArrowForceKernel::execute() line 306  pFrc[0].x = %f ", FORCE_TO_DBL(pFrc[0].x));
    // printf("pFrc[0].y = %f ", FORCE_TO_DBL(pFrc[0].y)); 
    // printf("pFrc[0].z = %f \n", FORCE_TO_DBL(pFrc[0].z));     

    double energy = 0;

    std::string platform_name = context.getPlatform().getName();
    int natoms = context.getSystem().getNumParticles();

    vector<Vec3> forces_loc(natoms);
    for (int i = 0; i < natoms; i++)
    {
        forces_loc[i][0] = FORCE_TO_DBL( pFrc[i].x ) * 41.84 * scale_force;
        forces_loc[i][1] = FORCE_TO_DBL( pFrc[i].y ) * 41.84 * scale_force;
        forces_loc[i][2] = FORCE_TO_DBL( pFrc[i].z ) * 41.84 * scale_force;
    }

    //for( int i = 0; i < 2; i++ )
    //{
    //    printf(" CommonCalcArrowForceKernel::execute() line 367  forces[%d] = %f %f %f \n", i, forces_loc[i][0], forces_loc[i][1], forces_loc[i][2]);
    //}
    //printf(" CommonCalcArrowForceKernel::execute() line 369 \n" );

    addForces(forces_loc, context);

    // Potential energy in kJ/mol in OpenMM
    energy = pEnv->m_pCmptMatrix->m_pfCalculatedValues[0][SimulationCore::ECalculatedValues::eEnergPot]*4.184*scale_force; 
    return energy;
}

void CommonCalcArrowForceKernel::copyParametersToContext(ContextImpl& context, const ArrowForce& force) {
    // printf(" CommonCalcArrowForceKernel::copyParametersToContext() \n");
    // fflush(stdout);

    SimulationCore::CSimulationEnvironment* pEnv = pSysLdr->GetSimulationenvironment();
    //bool bRes = copyCrdFromArbalestToContext(pEnv, context); // Copy coordinates from Arbalest environment to OpenMM context as coordinates may change to moving atoms to the central unit cell 
}

bool CommonCalcArrowForceKernel::posInternalChanged(ContextImpl& context) const
{
    printf(" CommonCalcArrowForceKernel::posInternalChanged() = %d \n", this->positions_changed);
    fflush(stdout);
    return this->positions_changed;
}

void CommonCalcArrowForceKernel::copyInternalPositionsToContext(ContextImpl& context)
{
    context.setPositions(positions_internal);
    // context.setPeriodicBoxVectors(a, b, c); No changing of periodic box vectors so far
}


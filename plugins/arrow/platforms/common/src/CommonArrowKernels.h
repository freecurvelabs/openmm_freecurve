#ifndef COMMON_ARROW_KERNELS_H_
#define COMMON_ARROW_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "openmm/ArrowKernels.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeArray.h"

// Arbalest includes:

#include "CommonLib/Typedefs.h"
#include "ParseCommandLine.h"
#include "SystemLoadingLib/SystemLoader.h"
#include "SimulationCoreLib/SimulationEnvironment.h"
#include "SimulationCoreLib/SimulationEnvironmentGpu.h"
#include "SimulationToolsLib/SimulationController.h"
#include "GpuLib/GpuManager.h"
#include "CommonLib/Timer.h"
#include "CommonLib/SubprocessRun.h"

namespace OpenMM {

/**
 * This kernel is invoked by ArrowForce to calculate the forces acting on the system and the energy of the system.
 */
class CommonCalcArrowForceKernel : public CalcArrowForceKernel {
public:
    CommonCalcArrowForceKernel(const std::string& name, const Platform& platform, ComputeContext* p_cc) :
            CalcArrowForceKernel(name, platform), p_cc(p_cc), kernel_platform(platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ArrowForce this kernel will be used for
     */
    void initialize(const System& system, const ArrowForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ArrowForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const ArrowForce& force);

    /**
      * Copy coordinates from OpenMM context to Arbalest environment.
      *
      * @param context       the context ( source )
      * @param pEnvReplica   Arbalest environment Replica ( destination )
    */
    bool copyCrdFromContextToArbalest(ContextImpl& context, SimulationCore::CEnvironmentReplica* pEnvReplica);

    /**
      * Copy coordinates from Arbalest environment to OpenMM context.
      *
      * @param pEnvReplica   Arbalest environment Replica ( source )
      * @param context       the context ( destination )
      * 
    */
    bool saveInternalPositions(SimulationCore::CEnvironmentReplica* pEnvReplica, ContextImpl& context );
  
    /**
      * check of coordinates in Arbalest environment changes in the last force call.
      * @param context       the context to choose the proper kernel
    */
    bool posInternalChanged(ContextImpl& context) const;

    /**
     * Copy internal positions to OpenMM context.
     *
     * @param context       the context to copy the positions to
     */
    void copyInternalPositionsToContext(ContextImpl& context);

    /**
     * Set forces to the context
     *
     * @param forces      the ArrowForce to copy the parameters from
     */
    void addForces(vector<Vec3>& forces, ContextImpl& context);

private:

  const Platform& kernel_platform;
  ComputeContext* p_cc; // Pointer to the ComputeContext
  ComputeArray particleParams;
  ComputeArray pairParams;
	
	double scale_force;  // scale force parameter
  bool positions_changed = false;
  std::vector<OpenMM::Vec3> positions_internal; // internal atom positions in Arbalest environment

	std::vector<int> particle;

  // Arbalest structures:
  SimulationTools::CSimulationController SimController;
  std::shared_ptr<SystemLoading::CSystemLoader> pSysLdr;
	std::shared_ptr<SimulationCore::CMDSchemeCommonOperations> spMDSchemeOperations;
  
};

} // namespace OpenMM

#endif /*COMMON_ARROW_KERNELS_H_*/

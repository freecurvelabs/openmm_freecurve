#ifndef REFERENCE_ARROW_KERNELS_H_
#define REFERENCE_ARROW_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "ReferencePlatform.h"
#include "openmm/ArrowKernels.h"
#include "openmm/Vec3.h"

 // Arbalest includes:

#include "ParseCommandLine.h"
#include "SystemLoadingLib/SystemLoader.h"
#include "SimulationCoreLib/SimulationEnvironment.h"
#include "SimulationCoreLib/SimulationEnvironmentGpu.h"
#include "SimulationToolsLib/SimulationController.h"
#include "GpuLib/GpuManager.h"
#include "CommonLib/Timer.h"
#include "CommonLib/SubprocessRun.h"

#include <utility>
#include <vector>

class ReferenceConstraintAlgorithm;

namespace OpenMM {

/**
 * This kernel is invoked by ArrowForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcArrowForceKernel : public CalcArrowForceKernel {
public:
    ReferenceCalcArrowForceKernel(const std::string& name, const Platform& platform) : CalcArrowForceKernel(name, platform) {
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

private:

};


} // namespace OpenMM

#endif /*REFERENCE_ARROW_KERNELS_H_*/

#ifndef ARROW_KERNELS_H_
#define ARROW_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package 
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "openmm/ArrowForce.h"
#include "openmm/Kernel.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This kernel is invoked by ArrowForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcArrowForceKernel : public KernelImpl {
public:
    static std::string Name() {
        return "CalcArrowForce";
    }
    CalcArrowForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the DrudeForce this kernel will be used for
     */
    virtual void initialize(const System& system, const ArrowForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the DrudeForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const ArrowForce& force) = 0;
};


} // namespace OpenMM

#endif /*ARROW_KERNELS_H_*/

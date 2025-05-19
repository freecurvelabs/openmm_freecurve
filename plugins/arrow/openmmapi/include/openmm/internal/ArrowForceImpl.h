#ifndef OPENMM_ARROWFORCEIMPL_H_
#define OPENMM_ARROWFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package 
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "openmm/ArrowForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of ArrowForce.
 */

class OPENMM_EXPORT_ARROW ArrowForceImpl : public ForceImpl {
public:
    ArrowForceImpl(const ArrowForce& owner);
    ~ArrowForceImpl();
    void initialize(ContextImpl& context);
    const ArrowForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
    std::vector<std::pair<int, int> > getBondedParticles() const;
    bool usesPeriodicBoundaryConditions() const { return true; }
private:
    const ArrowForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_ARROWFORCEIMPL_H_*/

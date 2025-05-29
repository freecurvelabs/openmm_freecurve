#ifndef OPENMM_ARROWFORCE_H_
#define OPENMM_ARROWFORCE_H_

/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include <iostream>
#include "internal/windowsExportArrow.h"

namespace OpenMM {

/**
 * This class implements forces that are specific to ARROW force field. 
 */

class OPENMM_EXPORT_ARROW ArrowForce : public Force {
public:
    /**
     * Create a ArrowForce.
     */
    ArrowForce();

    /**
     * Set name of the Arbalest config file to setup Arbalest structures
     */
    void setArbalestConfig(const std::string arb_config_fname) {
        this->arb_config_fname = arb_config_fname;
    }

    /**
     * Get name of the Arbalest config file to setup Arbalest structures
     */
    std::string getArbalestConfig() const {
        return arb_config_fname;
    }

    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const {
        return particles.size();
    }

    int addParticle(int particle);
    
    /**
     * Get the parameters for a Arrow particle.
     *
     * @param index                the index of the Arrow particle for which to get parameters
     * @param[out] particle        the index within the System of the Arrow particle
     */
    void getParticleParameters(int index, int& particle) const;
    
    /**
     * Set the parameters for a Arrow particle.
     *
     * @param index           the index of the Arrow particle for which to set parameters
     * @param particle        the index within the System of the Arrow particle
     */
    void setParticleParameters(int index, int particle);

    /**
     * Check if Internal Positions were changed on the last force call .
     *
     */
    bool posInternalChanged(Context& context);

    /**
     * Copy internal positions to OpenMM context.
     *
     * @param context       the context to copy the positions to
     */
    void copyInternalPositionsToContext(Context& context);

    /**
     * Update parameters in context
     *
     * @param context        the context 
     */
    void updateParametersInContext(Context& context);
    
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if nonbondedMethod uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
         return true;
    }

    /**
     * Set Force scale factor - to use in mixed hamiltonians
     *
     * @param scale        scale factor for the force 
     */
    void setScaleForce(double scale) { this->scale_force = scale; } 

protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    std::vector<ParticleInfo> particles;
    std::string arb_config_fname;
public:
    double scale_force; // Force scale factor - to build mixed hamiltonians 
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class ArrowForce::ParticleInfo {
public:
    int particle;
    ParticleInfo() {
        particle = -1;
    }
    ParticleInfo(int particle) :
        particle(particle) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_ARROWFORCE_H_*/

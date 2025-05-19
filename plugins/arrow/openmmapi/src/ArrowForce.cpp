/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/ArrowForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ArrowForceImpl.h"

using namespace OpenMM;
using namespace std;

ArrowForce::ArrowForce() {
    scale_force = 1.0;
}

int ArrowForce::addParticle(int particle) {
    particles.push_back(ParticleInfo(particle));
    return particles.size()-1;
}

void ArrowForce::getParticleParameters(int index, int& particle) const {
    ASSERT_VALID_INDEX(index, particles);
    particle = particles[index].particle;
}

void ArrowForce::setParticleParameters(int index, int particle) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].particle = particle;

}

ForceImpl* ArrowForce::createImpl() const {
    return new ArrowForceImpl(*this);
}

void ArrowForce::updateParametersInContext(Context& context) {
    dynamic_cast<ArrowForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

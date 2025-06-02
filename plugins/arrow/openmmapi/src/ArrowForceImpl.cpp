/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ArrowForceImpl.h"
#include "openmm/ArrowKernels.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>

using namespace OpenMM;
using namespace std;

ArrowForceImpl::ArrowForceImpl(const ArrowForce& owner) : owner(owner) {
}

ArrowForceImpl::~ArrowForceImpl() {
}

void ArrowForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcArrowForceKernel::Name(), context);
    const System& system = context.getSystem();

    // Check for errors in the specification of particles.

    set<int> usedParticles;
    for (int i = 0; i < owner.getNumParticles(); i++) {
        int particle[5];
        double charge, k, k2, k3;
        owner.getParticleParameters(i, particle[0]);
        for (int i = 0; i < 2; i++) {
            if (particle[i] < 0 || particle[i] >= system.getNumParticles()) {
                stringstream msg;
                msg << "ArrowForce: Illegal particle index: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
            if (usedParticles.find(particle[i]) != usedParticles.end()) {
                stringstream msg;
                msg << "ArrowForce: Particle index is used by two different Drude particles: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
            usedParticles.insert(particle[i]);
        }
        for (int i = 2; i < 5; i++) {
            if (particle[i] < -1 || particle[i] >= system.getNumParticles()) {
                stringstream msg;
                msg << "ArrowForce: Illegal particle index: ";
                msg << particle[i];
                throw OpenMMException(msg.str());
            }
        }
    }

    kernel.getAs<CalcArrowForceKernel>().initialize(context.getSystem(), owner);
}

double ArrowForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0) 
        return kernel.getAs<CalcArrowForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> ArrowForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcArrowForceKernel::Name());
    return names;
}

void ArrowForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcArrowForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

bool ArrowForceImpl::posInternalChanged(ContextImpl& context) {
    bool res = kernel.getAs<CalcArrowForceKernel>().posInternalChanged(context);
    return res;
}

OpenMM::Vec3 ArrowForceImpl::getArbalestShiftVec() const {
    OpenMM::Vec3 shift = kernel.getAs<CalcArrowForceKernel>().getArbalestShiftVec();
    return shift;
}

void ArrowForceImpl::copyInternalPositionsToContext(ContextImpl& context)
{
    kernel.getAs<CalcArrowForceKernel>().copyInternalPositionsToContext(context);
}

vector<pair<int, int> > ArrowForceImpl::getBondedParticles() const {
    int numParticles = owner.getNumParticles();
    vector<pair<int, int> > bonds(numParticles);
    for (int i = 0; i < numParticles; i++) {
        int p2, p3, p4;
        double charge, polarizability, aniso12, aniso34;
        owner.getParticleParameters(i, bonds[i].first);
    }
    return bonds;
}

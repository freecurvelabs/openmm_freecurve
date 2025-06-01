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


void ReferenceCalcArrowForceKernel::initialize(const System& system, const ArrowForce& force) {
    
}


double ReferenceCalcArrowForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    //cout << "ReferenceCalcArrowForceKernel::execute()" << std::endl;
    return 0.0; 
}

void ReferenceCalcArrowForceKernel::copyParametersToContext(ContextImpl& context, const ArrowForce& force) {
    
}

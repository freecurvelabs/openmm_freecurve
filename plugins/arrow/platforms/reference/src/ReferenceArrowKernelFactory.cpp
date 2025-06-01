/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "ReferenceArrowKernelFactory.h"
#include "CommonArrowKernels.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
            ReferenceArrowKernelFactory* factory = new ReferenceArrowKernelFactory();
            platform.registerKernelFactory(CalcArrowForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerArrowReferenceKernelFactories() {
    registerKernelFactories();
}

KernelImpl* ReferenceArrowKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    if (name == CalcArrowForceKernel::Name())
        return new CommonCalcArrowForceKernel(name, platform, NULL);

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}

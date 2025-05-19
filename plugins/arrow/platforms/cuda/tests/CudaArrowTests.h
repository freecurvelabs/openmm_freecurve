/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "CudaTests.h"

extern "C" void registerArrowCudaKernelFactories();

using namespace OpenMM;

void setupKernels (int argc, char* argv[]) {
    registerArrowCudaKernelFactories();
    platform = dynamic_cast<CudaPlatform&>(Platform::getPlatformByName("CUDA"));
    initializeTests(argc, argv);
}

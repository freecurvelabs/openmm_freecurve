/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "ReferenceTests.h"

extern "C" void registerArrowReferenceKernelFactories();

using namespace OpenMM;

void setupKernels (int argc, char* argv[]) {
    registerArrowReferenceKernelFactories();
    platform = dynamic_cast<ReferencePlatform&>(Platform::getPlatformByName("Reference"));
    initializeTests(argc, argv);
}

/* -------------------------------------------------------------------------- *
 *                                OpenMMArrow                                 *
 * -------------------------------------------------------------------------- */

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/ArrowForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace std;

//extern "C" void registerArrowSerializationProxies();

void testSerialization() {
    // Create a Force.

    // ArrowForce force1;
}

int main() {
    try {
//        registerArrowSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}


/* -------------------------------------------------------------------------- *
 *                               OpenMMArrow                                  *
 * -------------------------------------------------------------------------- *
 *    This a ARROW force field plugin to OpenMM package
 *    based on Arbalest program of Freecurve Inc
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/ArrowForce.h"
#include "SimTKOpenMMUtilities.h"
#include <iostream>
#include <vector>
#include "boost/format.hpp"

using namespace OpenMM;
using namespace std;

void testWaterDimer() {
    System system;
    system.addParticle(16.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(16.0);
    system.addParticle(1.0);
    system.addParticle(1.0);

    ArrowForce* arrow = new ArrowForce();

    std::cout << " testWaterDimer() pt 1 " << std::endl;
    //arrow->setArbalestConfig("wat_2_arrow_rerun_conf_templ_oo_R5.xml");
    arrow->setArbalestConfig("wat_2_arrow_calc_ene.xml");
    std::cout << " testWaterDimer() pt 2 " << std::endl;

    system.addForce(arrow);
    std::cout << " testWaterDimer() pt 3 " << std::endl;

    vector<Vec3> positions(6); // in nm
    Vec3 a, b, c;              // in nm
    positions[0] = Vec3( 0.000, 0.000, 0.000) * 0.1;
    positions[1] = Vec3( 0.957, 0.000, 0.000) * 0.1;
    positions[2] = Vec3(-0.240, 0.927, 0.000) * 0.1;
    positions[3] = Vec3( 3.140,-2.290, -0.960) * 0.1;
    positions[4] = Vec3( 3.770,-2.360,-0.280) * 0.1;
    positions[5] = Vec3( 2.780,-3.240,-1.160) * 0.1;

    a = Vec3( 3.2,  0.0,  0.0);  
    b = Vec3(   0,  3.2,  0.0);
    c = Vec3(   0,  0.0,  3.2);

    VerletIntegrator integ(1.0);

    Context context(system, integ, platform);
    std::cout << " testWaterDimer() pt 4 " << std::endl;

    context.setPositions(positions);
    context.setPeriodicBoxVectors(a, b, c);

    std::cout << " testWaterDimer() pt 5 " << std::endl;
    State state = context.getState(State::Energy | State::Forces);
    std::cout << " testWaterDimer() pt 6 " << std::endl;

    double ene_pot_arrow = state.getPotentialEnergy();
    std::cout << " testWaterDimer() pt 7 " << std::endl;
    std::vector<OpenMM::Vec3> forces = state.getForces();

    std::cout << "Computed Arrow Energy: " << ene_pot_arrow << std::endl;
    std::cout << "Computed Arrow Forces: " << std::endl;
    for (int i = 0; i < forces.size(); i++)
    {
        std::cout << boost::format("%3d %10.5f %10.5f %10.5f") % i % forces[i][0] % forces[i][1] % forces[i][2] << std::endl ;
    }
}

void setupKernels(int argc, char* argv[]);
void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        setupKernels(argc, argv);
        testWaterDimer();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "FAIL - ERROR.  Test failed." << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}

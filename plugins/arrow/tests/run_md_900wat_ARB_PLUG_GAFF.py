from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile("wat900.gro")
positions = gro.positions
a = openmm.Vec3(32.0,  0.0,  0.0) * 0.1
b = openmm.Vec3( 0.0, 32.0,  0.0) * 0.1
c = openmm.Vec3( 0.0,  0.0, 32.0) * 0.1

top = GromacsTopFile('900H2O_gmx.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
#system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=None, rigidWater=False)

system.removeForce(1)

frc_a = openmm.ArrowForce()
frc_a.setArbalestConfig("wat900_gaff_ene_arb.xml")
system.addForce(frc_a)

system.removeForce(0)
system.removeForce(0)

num_forces = system.getNumForces()
print(f"Number of forces in the system: {num_forces}")
for i in range(num_forces):
    force = system.getForce(i)
    print(f"Force {i}: {type(force).__name__}")
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
#integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.00025*picoseconds)
platform = openmm.Platform.getPlatformByName("Reference")
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)
simulation.context.setPeriodicBoxVectors(a, b, c)
state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print(f"Potential Energy: {potential_energy}")
#simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('wat900_arb_plugin_gaff_1fs_nohv.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, time=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(StateDataReporter("wat_md_arb_plugin_gaff_1fs_nohr.out", 10, step=True, time=True, potentialEnergy=True, temperature=True, density=True))
simulation.step(1000000)

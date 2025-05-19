from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile("wat900.gro")
positions = gro.positions

top = GromacsTopFile('900H2O_gmx.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=None,rigidWater=False)
#system.setDefaultPeriodicBoxVectors(a, b, c)

#system.removeForce(1)
system.removeForce(0)

frc_a = openmm.ArrowForce()
frc_a.setArbalestConfig("wat900_arb.xml")
system.addForce(frc_a)

#system.removeForce(0)
#system.removeForce(0)
#system.removeForce(0)

barostat = MonteCarloBarostat(1.0*atmospheres, 300.0*kelvin, 25)
system.addForce(barostat)

num_forces = system.getNumForces()
print(f"Number of forces in the system: {num_forces}")
for i in range(num_forces):
    force = system.getForce(i)
    print(f"Force {i}: {type(force).__name__} usesPeriodicBoundaryConditions() = {force.usesPeriodicBoundaryConditions()}")

#integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.0005*picoseconds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.0005*picoseconds)
platform = openmm.Platform.getPlatformByName("Reference")
#platform = openmm.Platform.getPlatformByName("Cuda")
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)
simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print(f"Potential Energy Init: {potential_energy}")
#simulation.minimizeEnergy()
state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print(f"Potential Energy Min: {potential_energy}")
simulation.reporters.append(DCDReporter('wat900_arrow_cp_05fs_cuda_rel.dcd', 20000))
simulation.reporters.append(StateDataReporter(stdout, 100, time=True,step=True, potentialEnergy=True, temperature=True,density=True))
simulation.reporters.append(StateDataReporter("wat900_md_cp_05fs_nohr_cuda_rel.out", 1000, step=True, potentialEnergy=True, temperature=True,density=True))
simulation.step(2000000)

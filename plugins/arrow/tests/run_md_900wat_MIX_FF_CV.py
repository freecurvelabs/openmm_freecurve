from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile("wat900_BOX30.gro")
positions = gro.positions

top = GromacsTopFile('900H2O_gmx.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=None,rigidWater=False)

#system.removeForce(1) # for rigidWater=False
#system.removeForce(0) # for rigidWater=True  

gaff_scale = 0.01 
gaff_nonbond_force = system.getForce(0)
for particle_index in range(gaff_nonbond_force.getNumParticles()):
    charge, sigma, epsilon = gaff_nonbond_force.getParticleParameters(particle_index)
    scaled_epsilon = epsilon * gaff_scale 
    scaled_charge = charge * gaff_scale 
    gaff_nonbond_force.setParticleParameters(particle_index, scaled_charge, sigma, scaled_epsilon)

frc_a = openmm.ArrowForce()
frc_a.setArbalestConfig("wat900_arb.xml")
frc_a.setScaleForce(0.99)
system.addForce(frc_a)

num_forces = system.getNumForces()
print(f"Number of forces in the system: {num_forces}")
for i in range(num_forces):
    force = system.getForce(i)
    print(f"Force {i}: {type(force).__name__} usesPeriodicBoundaryConditions() = {force.usesPeriodicBoundaryConditions()}")

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
#platform = openmm.Platform.getPlatformByName("Reference")
platform = openmm.Platform.getPlatformByName("CPU")
#platform.setPropertyDefaultValue('CpuThreads', '8')
simulation = Simulation(top.topology, system, integrator, platform)
#platform.setPropertyValue(simulation.context, "Threads", "8")
gaff_nonbond_force.updateParametersInContext(simulation.context)
simulation.context.setPositions(gro.positions)
simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
#simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('wat900_mixff_01_cv_2fs_fixedw_test.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 50, time=True,step=True, potentialEnergy=True, temperature=True,density=True))
simulation.reporters.append(StateDataReporter("wat900_mixff_01_cv_2fs_fixedw_test.out", 1000, step=True, potentialEnergy=True, temperature=True,density=True))
simulation.step(500000)

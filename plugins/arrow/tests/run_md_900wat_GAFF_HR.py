from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile('wat900_BOX30.gro')
top = GromacsTopFile('900H2O_gmx.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=None)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
#simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('wat900_openmm_gaff_1fs_hr_cv_test.dcd', 10000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, time=True,potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(StateDataReporter("wat_gaff_md_1fs_hr_cv_test.out", 100, time=True,step=True, potentialEnergy=True, temperature=True, density=True))
simulation.step(1000000)

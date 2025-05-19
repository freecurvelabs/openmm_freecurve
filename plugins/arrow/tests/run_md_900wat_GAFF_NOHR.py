from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile('900H2O.gro')
top = GromacsTopFile('900H2O_flex_gmx.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
#system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=None, rigidWater=False)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.00025*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
#simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('wat900_openmm_gaff_025fs_nohr.dcd', 4000))
simulation.reporters.append(StateDataReporter(stdout, 400, step=True, time=True,potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(StateDataReporter("wat_gaff_md_025fs_nohr.out", 40, time=True,step=True, potentialEnergy=True, temperature=True, density=True))
simulation.step(4000000)

# %%
import os
import openmm
import openmm.app
#os.chdir("C:\\MYPROG\\OPENMM_ARROW\\plugins\\arrow\\tests")
print(os.getcwd())
print("Loaded plugins:")
for plugin in openmm.pluginLoadedLibNames:
  print(plugin)

# %%
system = openmm.System()
for i in range(3):
  system.addParticle(16.0)
  system.addParticle(1.0)
  system.addParticle(1.0)
print(f"Natoms = {system.getNumParticles()}")

# %%
frc_a = openmm.ArrowForce() 
frc_a.setArbalestConfig("wat3_ene_arb.xml")
system.addForce(frc_a)

# %%
positions = []
positions.append(openmm.Vec3( 1.95 , 1.5 , 1.8 ) * 0.1)
positions.append(openmm.Vec3( 1.82 , 1.5 , 0.84) * 0.1)
positions.append(openmm.Vec3( 1.05 , 1.5 , 2.16) * 0.1)
positions.append(openmm.Vec3( 4.95 , 1.5 , 1.8) * 0.1)
positions.append(openmm.Vec3( 4.82 , 1.5 , 0.84) * 0.1)
positions.append(openmm.Vec3( 4.05 , 1.5 , 2.16) * 0.1)
positions.append(openmm.Vec3( 7.95 , 1.5 , 1.8) * 0.1)
positions.append(openmm.Vec3( 7.82 , 1.5 , 0.84) * 0.1)
positions.append(openmm.Vec3( 7.05 , 1.5 , 2.16) * 0.1)

a = openmm.Vec3(32.0,  0.0,  0.0) * 0.1
b = openmm.Vec3( 0.0, 32.0,  0.0) * 0.1
c = openmm.Vec3( 0.0,  0.0, 32.0) * 0.1


# %%
integrator = openmm.VerletIntegrator(1.0)
platform = openmm.Platform.getPlatformByName("Reference")
#platform.loadPluginsFromDirectory("C:\\CONDA_WIN_OPENMM_8.1.1_DEV\\lib")
#platform.loadPluginsFromDirectory("C:\\CONDA_WIN_OPENMM_8.1.1_DEV\\lib\\plugins")
#platform.loadPluginsFromDirectory("/home/igor/PROG_SRC/OPENMM_PYTHON_LINUX/lib")
#platform.loadPluginsFromDirectory("/home/igor/PROG_SRC/OPENMM_PYTHON_LINUX/lib/plugins")

context = openmm.Context(system, integrator, platform)
context.setPositions(positions)
context.setPeriodicBoxVectors(a, b, c)

state = context.getState(getEnergy=True,getForces=True)

ene=state.getPotentialEnergy()
frc=state.getForces()
print(f"\n Energy = {ene} \n")
print(f"Forces: ")
for i in range(len(frc)):
    print(f"{i:5d} {frc[i].x:12.6f} {frc[i].x:12.6f} {frc[i].x:12.6f} ")
    
#pdb = openmm.app.PDBFile('wat_2.pdb')
#simulation = openmm.app.Simulation(pdb.topology, system, integrator)
#simulation.context.setPositions(positions)
#simulation.context.SetPlatform(platform)
#simulation.minimizeEnergy()




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
nmol = 900
for i in range(nmol):
  system.addParticle(16.0)
  system.addParticle(1.0)
  system.addParticle(1.0)
print(f" Num Atoms = {system.getNumParticles()}")

# %%
frc_a = openmm.ArrowForce() 
frc_a.setArbalestConfig("wat900_arb.xml")
system.addForce(frc_a)

# %%
gro = openmm.app.GromacsGroFile("wat900.gro")
positions = gro.positions

a = openmm.Vec3(32.0,  0.0,  0.0)
b = openmm.Vec3( 0.0, 32.0,  0.0)
c = openmm.Vec3( 0.0,  0.0, 32.0)
#a = openmm.Vec3( 3.2,  0.0,  0.0)
#b = openmm.Vec3( 0.0,  3.2,  0.0)
#c = openmm.Vec3( 0.0,  0.0,  3.2)


# %%
integrator = openmm.VerletIntegrator(1.0)
platform = openmm.Platform.getPlatformByName("Reference")

context = openmm.Context(system, integrator, platform)
context.setPositions(positions)
context.setPeriodicBoxVectors(a, b, c)

state = context.getState(getEnergy=True,getForces=True)

ene=state.getPotentialEnergy()
frc=state.getForces()
print(f"\n Energy = {ene} \n")
print(f"Forces: ")
for i in range(5):
    print(f"{i:5d} {frc[i].x:12.3f} {frc[i].y:12.3f} {frc[i].z:12.3f} ")
    




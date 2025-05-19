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
system.addParticle(16.0)
system.addParticle(1.0)
system.addParticle(1.0)
system.addParticle(16.0)
system.addParticle(1.0)
system.addParticle(1.0)
system.getNumParticles()

# %%
frc_a = openmm.ArrowForce() 
frc_a.setArbalestConfig("wat_2_arrow_calc_ene.xml")
system.addForce(frc_a)

# %%
positions = []
positions.append(openmm.Vec3( 0.000, 0.000, 0.000) * 0.1)
positions.append(openmm.Vec3( 0.957, 0.000, 0.000) * 0.1)
positions.append(openmm.Vec3(-0.240, 0.927, 0.000) * 0.1)
positions.append(openmm.Vec3( 3.140,-2.290,-0.960) * 0.1)
positions.append(openmm.Vec3( 3.770,-2.360,-0.280) * 0.1)
positions.append(openmm.Vec3( 2.780,-3.240,-1.160) * 0.1)

a = openmm.Vec3(32.0,  0.0,  0.0) * 0.1
b = openmm.Vec3( 0.0, 32.0,  0.0) * 0.1
c = openmm.Vec3( 0.0,  0.0, 32.0) * 0.1


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
print(f"\n Energy in kcal/mol = {ene.value_in_unit(openmm.unit.kilocalories_per_mole)} \n")
print(f"Forces: ")
for i in range(len(frc)):
    print(f"{i:5d} {frc[i].x:12.6f} {frc[i].x:12.6f} {frc[i].x:12.6f} ")
    
#pdb = openmm.app.PDBFile('wat_2.pdb')
#simulation = openmm.app.Simulation(pdb.topology, system, integrator)
#simulation.context.setPositions(positions)
#simulation.context.SetPlatform(platform)
#simulation.minimizeEnergy()




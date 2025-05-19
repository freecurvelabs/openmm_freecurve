# %%
import os
import openmm
import openmm.app

# %%
system = openmm.System()
system.addParticle(16.0)
system.addParticle(1.0)
system.addParticle(1.0)

# %%
frc_a = openmm.ArrowForce() 
frc_a.setArbalestConfig("wat1_ene_arb.xml")
system.addForce(frc_a)

# %%
positions = []
positions.append(openmm.Vec3(  0.1000000000, -0.002337688334, 0.000) * 0.1)
positions.append(openmm.Vec3(  0.9571370208,  0.001729024256, 0.000) * 0.1)
positions.append(openmm.Vec3( -0.2379653731,  0.9273785949  , 0.000) * 0.1)

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
    print(f"{i:5d} {frc[i].x:12.6f} {frc[i].y:12.6f} {frc[i].z:12.6f} ")

ene0 = ene;
dx = 0.0001
dv = []
dv.append(openmm.Vec3(dx,  0.0,  0.0))  
dv.append(openmm.Vec3(0.0,  dx,  0.0))  
dv.append(openmm.Vec3(0.0,  0.0, dx ))  
for i in range(3):
  fc = [0.0,0.0,0.0]
  for j in range(3):
    positions[i] += dv[j]  
    context.setPositions(positions) 
    state = context.getState(getEnergy=True)
    ene1=state.getPotentialEnergy()
    positions[i] -= dv[j]  
    positions[i] -= dv[j]  
    context.setPositions(positions) 
    state = context.getState(getEnergy=True)
    ene2=state.getPotentialEnergy()
    dene = ene1 - ene2
    fc[j] = dene/(2.0*dx) 
    positions[i] += dv[j]  
  print( f" fc[{i}] = {fc[0]} {fc[1]} {fc[2]} ")





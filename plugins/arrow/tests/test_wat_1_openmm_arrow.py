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
#positions.append(openmm.Vec3( -0.0018852445, -0.002337688334, 0.000) * 0.1)
#positions.append(openmm.Vec3(  0.9571370208,  0.001729024256, 0.000) * 0.1)
#positions.append(openmm.Vec3( -0.2379653731,  0.9273785949  , 0.000) * 0.1)
# shifted as in Arbalest:
positions.append(openmm.Vec3(  15.63852893,  15.53514186, 16.0) * 0.1)
positions.append(openmm.Vec3(  16.5975512,   15.53920857, 16.0) * 0.1)
positions.append(openmm.Vec3(  15.4024488,   16.46485814, 16.0) * 0.1)

a = openmm.Vec3(32.0,  0.0,  0.0) * 0.1
b = openmm.Vec3( 0.0, 32.0,  0.0) * 0.1
c = openmm.Vec3( 0.0,  0.0, 32.0) * 0.1

# %%
integrator = openmm.VerletIntegrator(1.0)
#platform = openmm.Platform.getPlatformByName("Reference")
#platform = openmm.Platform.getPlatformByName("CUDA")
platform = openmm.Platform.getPlatformByName("CPU")

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





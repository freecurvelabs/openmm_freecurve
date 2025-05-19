from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

pdb = PDBFile('wat900.pdb')
positions = pdb.positions

print(positions[0])


# Print the positions of each atom
#for i, pos in enumerate(positions):
#    print(f"Atom {i+1}: {pos}")

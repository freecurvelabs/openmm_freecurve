#!/bin/bash 
gmx grompp -f md.mdp -c 900H2O.gro -p 900H2O_gmx.top -o wat900_md.tpr -maxwarn 30

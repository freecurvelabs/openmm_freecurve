#!/bin/bash 
gmx grompp -f md_wat_flex.mdp -c 900H2O.gro -p 900H2O_gmx.top -o wat900_md_flex.tpr -maxwarn 30

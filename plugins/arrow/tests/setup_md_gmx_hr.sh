#!/bin/bash 
gmx grompp -f md_hr.mdp -c 900H2O.gro -p 900H2O_gmx.top -o wat900_md_hr.tpr -maxwarn 30

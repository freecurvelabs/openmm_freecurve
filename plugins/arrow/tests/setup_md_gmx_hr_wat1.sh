#!/bin/bash 
gmx grompp -f md_hr.mdp -c wat_1.gro -p wat1_gmx.top -o wat1_md_hr.tpr -maxwarn 30

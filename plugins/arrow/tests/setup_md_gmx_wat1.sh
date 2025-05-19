#!/bin/bash 
gmx grompp -f md.mdp -c wat_1.gro -p wat1_gmx.top -o wat1_md.tpr -maxwarn 30

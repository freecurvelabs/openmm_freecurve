#!/bin/bash 
gmx mdrun -nt 4 -gpu_id 0 -s wat900_md_flex.tpr -deffnm WAT900_MD_FLEX     

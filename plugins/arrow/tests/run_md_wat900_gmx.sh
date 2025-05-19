#!/bin/bash 
gmx mdrun -nt 4 -gpu_id 0 -s wat900_md.tpr -deffnm WAT900_MD     

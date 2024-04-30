#!/bin/bash

#SBATCH -p parallel
#SBATCH -C skylake
#SBATCH -n 1                    
#SBATCH -c 64
#SBATCH -t 05:30:00              # Run time (d-hh:mm:ss)

##SBATCH -A m2_jgu-smrsim
#SBATCH -A m2_zdvhpc

srun hpc_inbreeding_script.R $@

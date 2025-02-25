#!/bin/tcsh

#BSUB -J strcutural_diversity		# Job name
#BSUB -o output.%J		            # Standard output file
#BSUB -e error.%J		            # Standard error file
#BSUB -W 48:00			            # Maximum time (hh:mm)
#BSUB -n 1			                # Number of MPI processes
#BSUB -R "span[hosts=1]"	        # Use n cores on 1 node
#BSUB -R "rusage[mem=80GB]"	        # Memory requirement (per core)
#BSUB -q cnr			            # The queue that I want it to join

# Activate conda env
conda activate /usr/local/usrapps/tcsi/LouisLANDIS/envs/louis_R_env

# Move to the work directory
cd /share/tcsi/lagoodal/R/Data/

# Rscript
Rscript Process_Trees_HPC.R

conda deactivate

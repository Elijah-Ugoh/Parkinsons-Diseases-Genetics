#!/bin/bash
#SBATCH --job-name=Mars2_update
#SBATCH --output=Mars2_update_%j.out
#SBATCH --error=Mars2_update_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=02:00:00
#SBATCH --partition=sens

#This job script runs the flip_ref_alt.py script and renders the output to the working dir

# write to standard output
cat $0

# Load thge required modules for python
ml GCC/13.3.0 SciPy-bundle/2024.05

# copy input files to temp dir
cp -p flip_ref_alt.py Score_Mars.txt ../../../MPBC_HRC_Rsq03_updated.bim $NAISS_TMP

#change to execution
cd $NAISS_TMP

# define the command
python flip_ref_alt.py

# rescue all outputs back to local working dir
cp -p unmatched_Mars_snps2.txt updated_Mars_scores2.txt $SLURM_SUBMIT_DIR 

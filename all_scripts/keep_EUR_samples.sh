#!/bin/bash
#SBATCH --job-name=keep_EUR_samples
#SBATCH --output=keep_EUR_samples_%j.out
#SBATCH --error=keep_EUR_samples_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=sens

#This job script uses the text file EUR_samples.txt to filter out all non-european samples from the KGP dataset in preparation for PRSice analysis
# NB: If using plink1.x, make to use a --keep file with both FID & IID.

# write plink output message to standard output
cat $0

# load plink module
module load GCC/11.3.0 OpenMPI/4.1.4
module load PLINK/2.00a3.7 

# copy input files to temp dir
cp -p KGP_autosomes_cor* EUR_samples.txt $NAISS_TMP

cd $NAISS_TMP

# define the plink command 
plink2 --bfile KGP_autosomes_cor --keep EUR_samples.txt --make-bed --out EUR_KGP_autosomes_cor --threads 25

# rescue output back to local working dir
cp -p EUR_KGP_autosomes_cor* $SLURM_SUBMIT_DIR 

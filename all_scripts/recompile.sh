#!/bin/bash
#SBATCH --job-name=recompile
#SBATCH --output=recompile_%j.out
#SBATCH --error=recompile_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=sens

#This job script recompiles the plinkbinary files using the --extract option and the NO_INDEL_EUR_KGP_autosomes_cor.txt file

# write to standard output
cat $0

# load plink module
module load GCC/11.3.0 OpenMPI/4.1.4
module load PLINK/2.00a3.7 

# copy input files to temp dir
cp -p EUR_KGP_autosomes_cor* NO_INDELS_EUR_KGP_autosomes_cor.txt $NAISS_TMP

cd $NAISS_TMP

# define the plink command 
plink --bfile EUR_KGP_autosomes_cor --extract NO_INDELS_EUR_KGP_autosomes_cor.txt --make-bed --out NO_INDELS_EUR_KGP_autosomes_cor --threads 25

# rescue output back to local working dir
cp -p NO_INDELS_EUR_KGP_autosomes_cor* $SLURM_SUBMIT_DIR 

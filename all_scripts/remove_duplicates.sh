#!/bin/bash
#SBATCH --job-name=remove_duplicates 
#SBATCH --output=remove_duplicates_%j.out
#SBATCH --error=remove_duplicates_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=sens

#job script for extracting out the duplicate SNP IDs and creating new plink binary files

# write to standard output
cat $0

# load plink module
module load GCC/11.3.0 OpenMPI/4.1.4
module load PLINK/2.00a3.7 

# copy input files to temp dir
cp -p autosomes_KGP.* multi_allelic_variants.txt $NAISS_TMP

cd $NAISS_TMP

# define the plink command 
plink --bfile autosomes_KGP --exclude multi_allelic_variants.txt --make-bed --out autosomes_KGP_final --threads 25 

# rescue output back to local working dir
cp -p autosomes_KGP_final.* $SLURM_SUBMIT_DIR 

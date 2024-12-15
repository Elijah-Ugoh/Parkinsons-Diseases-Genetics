#!/bin/bash
#SBATCH --job-name=update_name
#SBATCH --output=update_name_%j.out
#SBATCH --error=update_name_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=sens

#This job script helps to rename variants in the KGP dataset to chr:bp:ref:alt to ensure consistency for further analysis

# write to standard output
cat $0

# load plink module
module load GCC/11.3.0 OpenMPI/4.1.4
module load PLINK/2.00a3.7 

# copy input files to temp dir
cp -p autosomes_KGP_final* KGP_SNPS_noalleles.txt $NAISS_TMP

cd $NAISS_TMP

# define the plink command 
plink --bfile autosomes_KGP_final --update-name KGP_SNPS_noalleles.txt --make-bed --out KGP_autosomes_cor --threads 15

# rescue output back to local working dir
cp -p KGP_autosomes_cor* $SLURM_SUBMIT_DIR 

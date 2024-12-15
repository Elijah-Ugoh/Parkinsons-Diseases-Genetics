#!/bin/bash
#SBATCH --job-name=bgen_to_plink 
#SBATCH --output=bgen_plink_%j.out
#SBATCH --error=bgen_plink_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=sens

#job script for converting the KGP files from bgen format to plink binary format

# write to standard output
cat $0

# load plink module
module load GCC/11.3.0 OpenMPI/4.1.4
module load PLINK/2.00a3.7 

# copy input files to temp dir
cp -p all_phase3.* $NAISS_TMP

cd $NAISS_TMP

# define the plink command 
plink2 --pfile all_phase3 vzs --max-alleles 2 --chr 1-22 --allow-extra-chr --make-bed --out autosomes_KGP --threads 25 

# rescue output back to local working dir
cp -p autosomes_KGP* $SLURM_SUBMIT_DIR 

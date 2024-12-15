#!/bin/bash
#SBATCH --job-name=prsice
#SBATCH --output=prsice_%j.out
#SBATCH --error=prsice_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=sens

#This job script runs the PRSice software and renders the output to the working dir

# write to standard output
cat $0

# load plink and PRSice modules
ml GCC/11.3.0 OpenMPI/4.1.4 PLINK/2.00a3.7 PRSice/v2.3.5

# copy input files to temp dir
cp -p PRSice.R prsice_run_07_12_2024.valid base_data.txt MPBC_HRC_Rsq03_updated.* NO_INDELS_EUR_KGP_autosomes_cor* -r PRSice_linux $NAISS_TMP

#change to execution
cd $NAISS_TMP

# define the command
Rscript PRSice.R --dir $NAISS_TMP \
 --prsice PRSice_linux \
 --base base_data.txt \
 --snp SNP --bp BP --chr CHR --A1 A1 --A2 A2 --stat BETA --pvalue P \
 --target MPBC_HRC_Rsq03_updated \
 --beta \
 --binary-target T \
 --extract prsice_run_07_12_2024.valid \
 --ld NO_INDELS_EUR_KGP_autosomes_cor --score avg --out prsice_run_08_12_2024 \
 --clump-kb 250kb \
 --clump-p 1.000000 \
 --clump-r2 0.100000 \
 --seed 782260214 \
 --thread 25

## rescue output back to local working dir
cp -p prsice_run_08_12_2024* $SLURM_SUBMIT_DIR 

#!/bin/bash
#SBATCH --job-name=prsice
#SBATCH --output=prsice_%j.out
#SBATCH --error=prsice_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --partition=sens

#This job script runs the PRSice software and renders the output to the working dir

# write to standard output
cat $0

# Load the required modules for R, Plink, and PRSice(not needed if PRSice executables is copied to the cluster as below)
ml GCC/11.3.0 OpenMPI/4.1.4 PLINK/2.00a3.7 PRSice/v2.3.5

# copy input files to temp dir
cp -p PRSice.R base_data.txt MPBC_HRC_Rsq03_updated.* NO_INDELS_EUR_KGP_autosomes_cor* -r PRSice_linux $NAISS_TMP

#change to execution
cd $NAISS_TMP

# define the command
BASE_CMD="Rscript PRSice.R --dir $NAISS_TMP \
 --prsice PRSice_linux \
 --base base_data.txt \
 --snp SNP --bp BP --chr CHR --A1 A1 --A2 A2 --stat BETA --pvalue P \
 --target MPBC_HRC_Rsq03_updated \
 --beta \
 --binary-target T \
 --ld NO_INDELS_EUR_KGP_autosomes_cor --score avg --out prsice_run_17_12_2024 \
 --clump-kb 250kb \
 --clump-p 1.000000 \
 --clump-r2 0.100000 \
 --seed 782260214 \
 --print-snp \
 --thread 25"

# also add the following options if necessary
# --cov NEW_COVS.txt
# --cov-col SEX,AGE,@PC[1-5]
# --perm 12000

# First run
echo "Running PRSice (First run)..."
eval $BASE_CMD

# Check if PRSice detected duplicates in the base data and have outputed a .valid file
# the the same output prefix to track the file
VALID_FILE="prsice_run_17_12_2024.valid"
if [[ -f "$VALID_FILE" ]]; then
    echo "Duplicate/ambigous SNPs detected/removed and .valid file created by PRSice. Re-running PRSice with the --extract option..."
    BASE_CMD="$BASE_CMD --extract $VALID_FILE"
    eval $BASE_CMD
else
    echo "No duplicate or ambigous SNPs detetecd by PRSice. No need for re-run..."
fi

# rescue all outputs back to local working dir
cp -p prsice_run_17_12_2024* $SLURM_SUBMIT_DIR 

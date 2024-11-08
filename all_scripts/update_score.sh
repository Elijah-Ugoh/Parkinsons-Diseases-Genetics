#!/bin/bash
#
#SBATCH -J dbSNP_filtered
#SBATCH --partition=sens
#SBATCH -o dbSNP_filtered_%j.out
#SBATCH -e dbSNP_filtered_%j.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

# write to standard output 
cat $0

# copy input dta to node local disk
cp -p -r CHECKSUMS GCF_000001405.25.gz* $NAISS_TMP

# change to the execution directory
cd $NAISS_TMP

# run  command
zgrep -w -f rsIDs.txt GCF_000001405.25.gz  > dbSNP_filtered.txt

# comment out the un-needed part to see if the result is different
# | awk '{print $3":"$4":"$5":"$6"\t"$7"\t"$8}' > updated_score_file.txt

# rescue results to the submission directory
cp -p dbSNP_filtered.txt $SLURM_SUBMIT_DIR


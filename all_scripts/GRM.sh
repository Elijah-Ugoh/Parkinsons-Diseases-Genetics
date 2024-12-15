#!/bin/bash
#
#SBATCH -J GCTA_GRM
#SBATCH --partition=sens
#SBATCH -o GRM_%j.out
#SBATCH -e GRM_%j.err
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# THis script runs the GCTA software for calculating genetic relatedness matrix among 1864 indiduals
# using about 17 million snps obtained from the plink binary input files

# write to standard output 
cat $0

# load gcta module
module load GCC/12.3.0
module load GCTA/1.94.1

# copy input dta to node local disk
cp -p pruned_data.bim pruned_data.bed pruned_data.fam gcta64 $NAISS_TMP

# change to the execution directory
cd $NAISS_TMP

# run gcta to compute genetic relatedness matrix
# partition the GRM into 3 parts for resource efficiency
gcta64 --bfile pruned_data --make-grm-part 3 1 --out GRM --autosome --maf 0.05 --thread-num 5
gcta64 --bfile pruned_data --make-grm-part 3 2 --out GRM --autosome --maf 0.05 --thread-num 5
gcta64 --bfile pruned_data --make-grm-part 3 3 --out GRM --autosome --maf 0.05 --thread-num 5 

# merge all parts together 
cat GRM.part_3_*.grm.id > GRM.grm.id
cat GRM.part_3_*.grm.bin > GRM.grm.bin
cat GRM.part_3_*.grm.N.bin > GRM.grm.N.bin

# run the relatedness threshold cutoff
gcta64 --grm-cutoff 0.125 --grm GRM --out GRM_0125 --make-grm --thread-num 5

# rescue results to the submission directory
cp -p GRM* $SLURM_SUBMIT_DIR

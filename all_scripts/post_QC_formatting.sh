#!/bin/bash
#SBATCH --job-name=perl
#SBATCH --output=perl_%j.out
#SBATCH --error=perl_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=00:30:00
#SBATCH --partition=sens-xl

# write to standard output
cat $0

# copy input files to temp dir
cp -p FILTERED.test.* HRC-1000G-check-bim-v4.2.pl filter_variants_maf_r2_updated.R HRC.r1-1.GRCh37.wgs.mac5.sites.tab $NAISS_TMP

#change to execution
cd $NAISS_TMP

perl HRC-1000G-check-bim-v4.2.pl -b FILTERED.test.bim -f FILTERED.test.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

## rescue output back to local working dir
cp -p *.txt Run-plink.sh $SLURM_SUBMIT_DIR 

# Polygenic Risk Score Calculation
## Using 90 Significant Risk Variants
With the covariate and beta score files created, the Polygenic Risk Scores for all individuals in the study can be calculated using plnk. 

- The ```--score``` function, which is commonly used to compute a weighted sum of allele dosages, based on the scores in the ```metaanalyis90.txt``` file is included.
- ```metaanalyis90.txt``` contains the PD-associated risk variants and their effect sizes.
- The binary PLINK files (.bed, .bim, and .fam) containing genotype data of 1864 individuas.
- The covariate file contains the individaul IDs, age, age at diagnosis (AAD), disease status, education, and sex information.

```bash
# first covert the SNP IDs in the beta score file to the corresponding SNP IDs in the .bim file
python update_score_ID.py MPBC_HRC_Rsq03_updated.bim metaanalysis90.txt meta_GRS_updated.txt

# run PRS computation
mkdir processed_analysis_files/04_PRS
plink --bfile MPBC_HRC_Rsq03_updated --score meta_GRS_updated.txt --out processed_analysis_files/04_PRS/PRS_PD_MPBC
```

Output:
- PRS_PD_MPBC.profile: Contains the GRS for each sample, calculated using valid variants.
- PRS_PD_MPBC.log: Log file
- For this analysis, all the variants from the score file matched those in the genotype data.

The SCORE column in the .profile file shows the computed risk score for each individual sample based on the SNPs present in the individual.

### Visualization in R
The PRS visualization is completed in R (R/4.3.2) using the ```PRS_PD_MPBC.profile``` and ```covariates.txt``` files and the ```PRS_analysis.R``` script.

```bash
# load the R module and its dependencies
moddule spider R/4.3.2
module load GCC/12.3.0
moddule load R/4.3.2

R < PRS_analysis.R --no-save
```

Separation noticed in the PRS for both the cases and controls was not strong.

The same analysis is also conducted for using the 1805 SNPs associated with Parkinson's disease in the Nall et al., 2019 study. 

## PRS Calcultion Uisng 1805 Risk Variants
Conidering that the SNPs file only has ```rsIDs```, this must be converted to the ```chrom:pos:ref:alt``` to align with the format in the .bim file before the analysis is done. 

```bash
# on the cluster, extract the rsIDs from the 1805 SNPs score file
mkdir processed_analysis_files/04_PRS/dbSNP/
cd processed_analysis_files/04_PRS/dbSNP/
cut -f 1 ../supplementary_data/1800snps_rs\ \(1\).tab > rsIDs.txt
```
- Here, rsIDs.txt is a file with only the rsIDs from the score file, ensuring we only pull matching entries from the dbSNP database.

NB: 
- Pulling directly from the dbSNP database did not yield the desired result, as it returned multiple alternate allele for each SNP. 
- The SNP ID conversion was done by comparing the rsIDs in the 1805 SNPs with the data in the complete Nalls Summary Statistics.

Next Steps:

- Download complete Nalls Summary Statistics with all SNPs data

```bash
# update summary statistics file name for easy handling
mv summary_statistics_Nalls_et_al_2019_EUR_metaGWAS_no23andme_hg38.txt updated_score_file.txt
```
Next, the file is re-formatted according to the structure in our .bim file
```bash 
# create a new file
echo -e "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight" > new_score_file.txt
# grab and apped re-arranged fields to the file
awk -F'[:| \t]' 'NR >1{
# Extract fields
chr = ($1 ~ /^chr/) ? substr($1, 4) : $1; # Remove "chr" to get the chromosome number and postion
pos = $2; # Extract position
rsID = $17; # Extract rsID
effect_size = $9; # Extract effect weight
effect_allele = toupper($3); # Extract effect allele
other_allele = toupper($4); # Extract other allele 
# Print in the required format
print rsID "\t" chr "\t" pos "\t" effect_allele "\t" other_allele "\t" effect_size
}' updated_score_file.txt >> new_score_file.txt

# update SNPs in score file to ensure they match the .bim data (SNP ID, effect, and alternate/other alleles)
sbatch update_scores_files.sh 

NB: The ```re_format_scores.py``` script is used for this. 

# Next, ensure the SNPs to be used are exactly the 1805 we need
awk 'NR==FNR {rsids[$1]; next} $4 in rsids { print $1, $2, $3}' rsIDs.txt updated_1805_scores.txt > final_score.txt
```

```bash
# PRS calculation can now be run again in PLINK
mkdir && cd rerun_1805/
module load GCC/11.3.0  OpenMPI/4.1.4 PLINK/2.00a3.7
plink --bfile ../../../../MPBC_HRC_Rsq03_updated --score ../final_score.txt --out rerun_1805
```

Output: 
- Only 1788 of the 1805 snps matched in the summary statistics
- All 1805 snps were loaded as valid predictors in PLINK

### Visualization in R
Again, the visualization is done in R (R/4.3.2), but this time, using the ```rerun_1805.profile``` file, and the ```covariates.txt``` (same as with the 90 variants) file.

```bash
# Modify the path to the variables (covariates.txt and 1800.profile) in the script as necessary
R < PRS_analysis.R --no-save
```  

## Using PRSice
To run PRS calculation using PRSice, the following files are required:
- Published GWAS or summary statistic or base data
- QC'ed, Imputed target data (plink files for the MPBC cohort)
- LD reference dataset (necessary if sample data is less than 500 individuals)
- Covarariates file (optional)

### Set Up
1. First download the relevant 1000 genome phase 3 database (build 37) from [here](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg) for LD referencing. 

```bash
wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1
zstd -d all_phase3.pgen.zst  # decompress the file
wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1
wget https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1
mv phase3_corrected.psam all_phase3.psam  # rename before use

# NB: wget may not work with dropbox, so download the files manually from a web browser. 
```

2. Next, convert the 1KGP data to plink binary format.

The plink command below converts the decompressed 1KGP data into PLINK binary format (.bed, .bim, .fam) for autosomal chromosomes (1–22). It excludes variants with more than two alleles (--max-alleles 2) and allows for extra chromosomes using ```--allow-extra-chr```.

```bash
# Due to meemory requiremnets, use the SLURM job script below to run the command on CSENS compute node 
# This applies to the rest of the commands executed leading up to and running PRSice.

```bash
# NB: In this scrpt, PLINK2 is used for handling BGEN files
sbatch bgen_plink.sh
```

3. Remove multi-allelic variants

After converting from BGEN to PLINK using the ```bgen_plink.sh``` slurm script, we extract and save the multi-allelic variants that have duplicated SNP IDs to a file.

```bash
# save a list of all SNPs with unique IDs
awk '{print $2}' autosomes_KGP.bim | sort | uniq -d > multi_allelic_variants.txt
```

This file is used along with the output from step 3 to produce new plink binaries: ```autosomes_KGP_final```
```bash
# use the script below
sbatch remove_duplictes.sh
```

4. Remove INDELS and Rename all variants to harmonize the dataset for downstream analyses

NB: The outputs from step 4 include INDELS and SNPS for all individuals
```bash
# so, we remove INDELS first, and leave only the SNPs
awk 'BEGIN {OFS="\t"} length($5) == 1 && length($6) == 1 {print}' autosomes_KGP_final.bim > no_indels_autosomes_KGP.txt

# next, update the IDs properly 
awk 'BEGIN {OFS="\t"} {print $2, $1":"$4":"$5":"$6}' no_indels_autosomes_KGP.txt > KGP_SNPS_noalleles.txt
```

```bash
# then rename all variants in the KGP dataset from rsIDs to the chr:bp:ref:alt format to match the MPBC data and GWAS summary stat
sbatch update_name.sh
```
This results in KGP_autosomes_cor, the final, corrected dataset with renamed variants in PLINK binary format. It serves as the 1000 Genomes autosomal reference panel for ancestry estimation and LD-based analyses.

5. Filter out only European samples in the KGP genotype data. This is beacuse we only want Europeans in the LD reference for PRSice.  

```bash
# separate out the European sample IDs in a different file
cat all_phase3.psam | grep 'EUR' | cut -f 1 > EUR_samples.txt

# with the output above, extract genotype data for only European samples
sbatch keep_EUR_samples.sh

# remove variants with INDELS or any multi-allelic representations from the output
awk 'BEGIN {OFS="\t"} length($5) == 1 && length($6) == 1 {print}' EUR_KGP_autosomes_cor.bim > NO_INDELS_EUR_KGP_autosomes_cor.txt
```
NB: It is necessary to retain only the SNPs because INDELS would create ambiguity and errors downstream when we run PRSice.

6. Recompile the plink binaries to ensure consistency with the data contained in the .bed and .fam files

```bash
# use the recompile.sh script. 
sbatch recompile.sh
```

7. Format the GWAS summary stat to include the fields needs needed for PRSice

Note that all SNPs are included and the summary stat has been formatted with these fields: 
```bash
CHR     BP      SNP     A1      A2      BETA        P      
```
- Previoulsy, the score files (90 and 1805 risk variants) had no p-values and only the effect alleles were included. 
- We still use the Nalls et al summary stats, wchich can be downloaded from the [PGS catalog](https://www.pgscatalog.org/publication/PGP000235/).

```bash
# move back to the 04_PRS dir
cd ../04_PRS
# extract the 'CHR' and 'BP' fields
awk 'NR>1 {print $1}' updated_score_file.txt | cut -f 2 -d ':' > ../../06_PRSice/bp.txt
awk 'NR>1 {print $1}' updated_score_file.txt | cut -f 1 -d ':' > ../../06_PRSice/chr.txt
sed -i 's/chr//g' ../../06_PRSice/chr.txt
paste  ../../06_PRSice/chr.txt  ../../06_PRSice/bp.txt > ../../06_PRSice/chr_bp.txt

# now back to the 06_PRSice folder
cd ../06_PRSice
# copy th raw summary stat to the cwd
cp ../dbSNP/updated_score_file.txt .
# extract relevant columns
awk 'NR>1 {print $1":"toupper($2)":"toupper($3)"\t"toupper($2)"\t"toupper($3)"\t"$8"\t"$10}' ../04_PRS/dbSNP/updated_score_file.txt > all_snps.txt
# remove all occurences of 'chr' in the SNP IDs
sed -i 's/chr//g' all_snps.txt
# create a header file
nano headers.txt 
# merge the file created earlier with the rest of the data
paste chr_bp.txt all_snps.txt > merged_score_file.txt
# add headers
cat headers.txt merged_score_file.txt > base_data.txt
```

NB: The 1000 Genome database was used for LD clumping in PRSice, but the PRS scores generated was only slighly better. The differnce is almost negligible, and this makes sense since the MPBC data has 1864 samples. 

### Run PRSice
A script, ```PRSice.sh``` is also used for this. First, create a new dir in the current dir:
- Then download and move the PRSice executable 'PRSice_linx' and the R script 'PRSice.R' to it. 
- Copy the base data, taget data, and QC'ed KPG data, ```NO_INDELS_EUR_KGP_autosomes_cor```  into the dir.
- Run PRSice as follows:
  
```bash
sbatch PRSice.sh

# PRSice automatically removes and saves duplicate SNPs from the summary stat, so the script is written to automatically re-run PRSice adding the --extract flag if duplicates are detected. The .valid contains the duplicates in this case.

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
```

To ensure, we have the best-fit regression without overfitting, we use empirical p-value by adding the ```--perm``` flag with 12000 permutation. Additionally, we can adjust for covariates in the regression using our covariates file by adding the ```--cov```and ```--cov-col``` options. See ```PRSice.sh``` script for full command.


Using the ```--perm``` option or adjusting for covariates didn't yield a better result. Moreover, this will be done in R when analysing the model fitness using logistic regression.

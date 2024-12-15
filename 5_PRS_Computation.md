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

No significant separation was noticed in the PRS for both the cases and controls.

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

# extract relevant columns from the SNPs data
awk 'NR>1 {print $1":"toupper($2)":"toupper($3)"\t"toupper($3)"\t"$8"\t"$16}' updated_score_file.txt > all_snps.txt

# remove 'chr' from the snp ids
sed -i 's/^chr//g' all_snps.txt

# next, extract the matching snps (1805) specified in the rsID.txt file above
awk 'NR==FNR { rsids[$1]; next} $4 in rsids { print $1, $2, $3}' rsIDs.txt all_snps.txt > cleaned_score_file.txt

# explore that final data
less cleaned_score_file.txt
wc -l cleaned_score_file.txt

# sort the chromosome numbers in ascending order (optional)
sort -t: -k1,1n -k2,2n cleaned_score_file.txt > cleaned_score_file2.txt
```

```bash
# run prs calculation in plink
plink --bfile ../MPBC_HRC_Rsq03_updated --score cleaned_score_file.txt --out 1805/1800
```

Alternatively, the steps above can also be run from the terminal using the python script ```update_IDs.py```
For this, a different version of the rsIDs file is needed. This will include the first 3 columns

```bash
cut -f 1-3 ../supplementary_data/1800snps_rs\ \(1\).tab > new_rsIDs.txt

# run the python script
python update_IDs.py

# next, run plink as above
```

Finally, using all the 7 million snps in summary statistics for the prs calculation

NB: plink terminates when runniung due to duplicate snp IDs in the ```all_snps.txt``` file. 
The python script below removes the duplicates

```python
import pandas as pd

# if pandas in not installed
conda install pandas 
pip install pandas

# Load the score file
score_df = pd.read_csv('filtered_snps.txt', sep='\t', header=None, names=[])

# Drop duplicate rows based on the 'chr_pos_ref_alt' column
score_df = score_df.drop_duplicates(subset=[0])

# Save the cleaned score file
score_df.to_csv('cleaned_score_file.txt', sep='\t', index=False, header=False)
```

```bash
# run prs calculation in plink
plink --bfile ../MPBC_HRC_Rsq03_updated --score cleaned_score_file.txt --out 1805/1800
```

Output: 
- Only 1788 of the 1805 snps matched in the summary statistics
- Of these, only 751 snps were valid predictors

### Visualization in R
Again, the visualization is done in R (R/4.3.2), but this time, using the ```1800.profile``` file, and the ```covariates.txt``` file.

```bash
# Modify the path to the variables (covariates.txt and 1800.profile) in the script as necessary
R < PRS_analysis.R --no-save
```
Similarly, there isn't significant difference in the cases and controls based on the 751 snp predicting PD.  

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

# PRSixe automatically removes and saves duplicate SNPs from the summary stat, so re-run it as shown below with the --extract flag if duplicates are detected

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



NB: Include QC for GWAS
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

## PRS Calcultion Uisng 1806 Risk Variants
conidering that the file only has ```rsIDs```, this must be converted to the ```chrom:pos:ref:alt``` to align with the format in the .bim file before the analysis is done. 

```bash
# first, a local version of the dbSNP database is created from "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/"

# make a new dir in the cosmos-home dir on local computer
mkdir dbSNP && cd dbSNP
# download dbSNP
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.md5
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi.md5
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/CHECKSUMS

# upload all to the working dir on CSENS
sftp elugoh@cs-diode.lunarc.lu.se
sftp> put -r dbSNP/ /home/elugoh/Documents/DataFolder/Elijah/

# on the cluster, extract the rsIDs from the 1805 SNPs score file
cd dbSNP/
cut -f 1 ../supplementary_data/1800snps_rs\ \(1\).tab > rsIDs.txt
```
- Here, rsIDs.txt is a file with only the rsIDs from the score file, ensuring we only pull matching entries from the dbSNP database.

- Next, run the ID sorting command as a job using optimal resources t reduce processing time.
```sbatch update_score.sh```

```bash
#!/bin/bash
#
#SBATCH -J update_score_file
#SBATCH --partition=sens
#SBATCH -o update_score_file_%j.out
#SBATCH -e update_score_file_%j.err
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

# write to standard output 
cat $0

# copy input dta to node local disk
cp -p -r CHECKSUMS GCF_000001405.25.gz.* $NAISS_TMP

# change to the execution directory
cd $NAISS_TMP

# run  command
zgrep -w -f rsIDs.txt GCF_000001405.25.gz  > dbSNP_filtered.txt

# comment out the un-needed part to see if the result is different
# | awk '{print $3":"$4":"$5":"$6"\t"$7"\t"$8}' > updated_score_file.txt

# rescue results to the submission directory
cp -p dbSNP_filtered.txt $SLURM_SUBMIT_DIR
```

- Process dbSNP entries in Python: Load the dbSNP_filtered.txt to create the mapping dictionary for the rsIDs, and replace each with the ```chrom:pos:ref:alt``` format if found.


NB: 
- Pulling directly from the dbSNP database didn not yield the desired result, as it returned multiple alternate allele for each SNP. 
- The SNP ID conversion was done by comparing the rsIDs in the 1805 SNPs with the data in the complete Nalls Summary Statistics.

Next Steps:

Download complete Nalls Summary Statistics with all SNPs data

```bash
# update summary statistics file name for easy handling
mv summary_statistics_Nalls_et_al_2019_EUR_metaGWAS_no23andme_hg38.txt updated_score_file.txt

# extract relevant columns from the data
awk 'NR>1 {print $1":"$2":"toupper($2)":"toupper($3)"\t"toupper($3)"\t"$8"\t"$16}' updated_score_file.txt > all_snps.txt

# remove 'chr' from the snp ids
sed -i 's/^chr//g' all_snps.txt

# next, extract the matching snps (1805) specified in the rsID.txt file above
awk 'NR==FNR { rsids[$1]; next} $4 in rsids { print $1, $2, $3}' rsIDs.txt all_snps.txt > fcleaned_score_file.txt

# explore that final data
less cleaned_score_file.txt
wc -l cleaned_score_file.txt

# sort the chromosome numbers in ascending order (optional)
sort -t: -k1,1n -k2,2n cleaned_score_file.txt > cleaned_score_file2.txt
```

```bash
# run prs calculation in plink
plink --bfile ../MPBC_HRC_Rsq03_updated --score cleaned_score_file.txt --out 1800/1800
```

Alternatively, this can also be run from the terminal using the python script ```update_IDs.py```
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
plink --bfile ../MPBC_HRC_Rsq03_updated --score cleaned_score_file.txt --out 1800/1800
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


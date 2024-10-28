# PRS Calculation
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

## Visualization in R
The PRS visualization is completed in R (R/4.3.2) using the ```PRS_PD_MPBC.profile``` and ```covariates.txt``` files and the ```PRS_analysis.R``` script.

```bash
R < PRS_analysis.R --no-save
```






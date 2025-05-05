# Creating a Covariate File

To run the polygenic risk score calculation and analysis, the following data files are needed:
- The QC-ed and Imputed genotype data files for all the individuals (cases and controls) in plink binary format (.bed, .bim & .fam files).
- A covariate file containing supplementary data about each individual, as well as their corresponding principal component analysis data.
- Finally, the score file containing the beta scores of each PD SNP of interest.

## Preparing the Covariate File
The information needed to prepare the covariate file is obtained from the questionaire issued to the individuals included in the MPBC cohort. Using a header description file, the appropriate columns were extracted and used to prepare the covariates.

NB: Not all the individuals in the cohort participated in the survey. 

### Extract Relevant Data from Questionnaire Data
```bash
# make new dir within Documents dir in the HPC to store project files
mkdir Elijah && cd Elijah
```

The ```QuestionnaireData_N1864_FINAL_CLEANED_210621.csv``` is in spreadsheet format, so it's converted to csv.

```bash
brew install gnumeric # to install xlsx-to-csv converter on mac
mkdir processed_analysis_files # create new dir within the project dir to store all new generated files
pwd # display cwd
"/home/elugoh/Documents/DataFolder/Elijah/"

head -1 QuestionnaireData_N1864_FINAL_CLEANED_210621_nomed.csv | tr "," "\n" | nl # view the file headers in new lines

# make another directory to hold the covariate file
cd processed_analysis_files && mkdir 02_covariate_file

# extract the columns/headers of interest from the questionnaire data to create a covariate file
cut -d ',' -f 2,4,6,10,179,181 QuestionnaireData_N1864_FINAL_CLEANED_210621_nomed.csv | sed 's/,/\t/g' | sed 's/"//g' > 02_covariate_file/covs.txt

# sed 's/,/\t/g' replaces commas with tabs
# sed -i "s/[\"-]//g' removes all the double quotes from the IDs and headers and the hyphens from the IDs

less -S 02_covariate_file/covs.txt # view new file
```

- The headers in the covariate file and the fam file are different and must be harmonized.

- This is done using an ID_match file that contains all IDs in the fam file matched to the corresponsing IDs in the questionnaire data/covariate file.

- We create an ID match file (first column matches the ID format in the fam file and second column matches the ID format in covariates file)

- The ID match file provided is saved with a trailing white space, which must be removed before it is used on the Linux terminal to ensure accurate ID matching results.

```bash
cd ..
cp ID_match.txt ID_match_unix.txt # make a copy 
dos2unix ID_match_unix.txt # convert the copy to unix format  
```

Now, create an updated covariate file using the corrected ID file.
```bash 
awk 'NR==FNR{a[$2]=$1; next} {if ($1 in a) $1=a[$1]; print}' ID_match_unix.txt processed_analysis_files/02_covariate_file/covs.txt > processed_analysis_files/02_covariate_file/new_cov.txt

# also add the familial ID (same as individual ID) to the file
cd processed_analysis_files/
awk 'BEGIN {OFS="\t"; print "FID"} NR>1 {print $1}' 02_covariate_file/new_cov.txt > 02_covariate_file/FID.txt # extract the FID column
paste 02_covariate_file/FID.txt 02_covariate_file/new_cov.txt > 02_covariate_file/updated_cov.txt # merge both files

# rename headers
cd 02_covariate_file
awk 'NR==1{$2="IID"; $4="SEX"; $5="EDUCATION"; $6="AAD"; $7="AGE"; print} NR>1' updated_cov.txt > renamed_updated_cov.txt
```

### Compute PCA
A principal component analysis is computed for all indiduals in the study to determine how many PCs to include in the covariates.

This is done on the filtered, pre-imputation data to get uninflated population stratification. 

To keep the working dir clean, create a different dir for the PCA calculation

```bash
cd ../ && mkdir 03_PCA
cd 03_PCA

plink --bfile ../01_QC/hwe/FILTERED.test --maf 0.05 --geno 0.01 --hwe 1E-6 --make-bed --out FILTERED.test_for_PCA

# Note that we still filter for genotyping, minor allele frequency, and hardy-weinberg equilibrium

plink --bfile FILTERED.test_for_PCA --indep-pairwise 50 5 0.5 --out prune 
plink --bfile FILTERED.test_for_PCA --extract prune.prune.in --make-bed --out prune
plink --bfile prune --pca --out NEW_MPBC_PCA
```
To determine the number of PCs to include in the covariates and GWAS analysis, we use the R scripts below to make a scree plot and a bar plot to see the variance explained by the eigenvalues of the PCs.

```bash
module load GCC/12.3.0 R/4.3.2
R < make_scree_plots.R --no-save
R < run_PCA.R --no-save
```

- ```make_scree_plots.R``` plot the variance explained by the eigenvalues in of all the PCs, as well as the first 10.  
- The R script, ```run_PCA.R``` can also be use to merge all the PCs and their individual IDs and saved as a text file in the current working dir as follows:
- The outputed table, ```pca_prunned_updated.txt``` has 20 PCs by default.
- Also, from the ```screePlot_MPBC_1-10.jpg``` and ```percentage_variance_explained.png``` plots, the first 10 PCs separate the data the most.

### Update Covariate File with PCs
Now, the covariate file can be updated with the right number of PCs needed. 
```bash
cd ../02_covariate_file/
# merge the covariates and the the first 10 PCs
awk 'NR==FNR {a[$1] = $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11; next} {print $0, (a[$1] ? a[$1] : "")}' ../03_PCA/pca_prunned_updated.txt renamed_updated_cov.txt > merged_pcs_cov.txt

# add the header names for the PCs
awk 'NR==1{$8="PC1"; $9="PC2"; $10="PC3"; $11="PC4"; $12="PC5"; $13="PC6";$14="PC7"; $15="PC8"; $16="PC9"; $17="PC10"; print} NR>1' merged_pcs_cov.txt > final_covariates.txt

# replace spaces in the documnet with tabs
sed -i 's/ /\t/g' final_covariates.txt

# Add maternal and paternal IDS to the covariate file to ensure the first 6 columns in both the covariate and QC'ed .fam files match exactly

cut -d " " -f 3-4 ../../MPBC_HRC_Rsq03_updated.fam > MAT_PAT_ID.txt # extract the columns from .fam file
awk 'BEGIN{print "MAT PAT"}1' MAT_PAT_ID.txt > MAT_PAT_ID2.txt # add headers
sed -i 's/ /\t/g' MAT_PAT_ID2.txt # replace spaces with tabs
paste final_covariates.txt MAT_PAT_ID2.txt > with_MAT_PAT_ID.txt # merge the new columns with the final covariate file

# re-arrange the columns
awk 'NR==FNR {cols[NR]=$18 "\t" $19; print $1 "\t" $2 "\t" cols[FNR] "\t" $4 "\t" $3 "\t" $5 "\t" $6"\t" $7 "\t" $8 "\t" $9"\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17}' c.txt > covariates.txt

less covariates.txt # view final coviarate file and compare the first two columns with the fam file
```

Finally, to ensure the IDs in the .fam file and covariate files match, we check for any differences using awk on the terminal
```bash
awk 'FNR==NR{a[$1];next} !($1 in a){print}' ../../MPBC_HRC_Rsq03_updated.fam covariates.txt
# If everything worked fine, the only difference should be the headers, which is missing in the .fam file.
```
```bash
FID     IID     MAT     PAT     SEX     PHENO       EDUCATION       AAD     AGE     PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10     
```


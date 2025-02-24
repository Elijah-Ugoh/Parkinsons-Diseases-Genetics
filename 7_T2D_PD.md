# Investigating the Impact of Type II Diabetes Mellitus (T2D) on Parkinson's Disease (PD) Risk in the MPBC Cohort
Steps:
- Data QC (extensive QC has been done on the MPBC cohort. So, no need for QC here, since the same dataset will be used in this pipeline).
- Create a new covariate file with additional information on T2D for all participants.
  - To create the covariate file, we must first extract diagnostic information from government data to ensure we include only T2D (and not T1D) in the analysis.
  - This is because clinical information in the questionnaire is self-reported and does not specify what type of diabetes.   
- Get summary stats from publicly available GWAS data and update them to plink-ready format.
- Match the SNPs in the summary stat with those in our QC'ed and imputed genotyping data to ensure alleles and ID consistency. 
- Run the PRS calculations using plink2. 
- Analyze the results using logistic regressisons and visualizations in R. 

## Creating Updated Covariate File with T2D Information
The steps below create an updated covariate file (using diabetes and parkinsons disease diagnoses from the Swedish national patient registry):

```bash
# create a wd in the pwd
mkdir 08_T2D_P && cd 08_T2D_PD

# move the diagnostic data to the cluster
>sftp put -r General_diagnostic_info/ /home/elugoh/Documents/DataFolder/Elijah/processed_analysis_files/08_T2D_PD/

# view the file to see which columns must be extracted from the questionnaire. 
cd General_diagnostic_info/ 
less -S UT_R_PAR_OV_14691_2021.txt | head -1 | tr "\t" "\n" | nl

# the individuals with either or both E11 (T2D) of G20 (PD) are filtered out to a separate file
(head -n 1 UT_R_PAR_OV_14691_2021.txt && grep -E 'E11|G20' UT_R_PAR_OV_14691_2021.txt) > G20_E11.txt # this keeps the headers
awk 'BEGIN {FS=OFS="\t"} {if ($1 == "") $1="NA"; print}' G20_E11.txt > G20_E11_NA.txt # assign NAs to first columns (main diagnosis) with missing entries
```
```bash
# print the relevant columns (IDs and diagnoses only) and redirect the output to a new file
awk 'BEGIN {FS=OFS="\t"}
NR ==1 {print $2, $3, "PD", "T2D"; next}
{
   has_G20 = 0;
   has_E11 = 0;
   # check first column (main diagnosis)
   if ($1 == "G20") has_G20 = 1;
   if ($1 == "E11") has_E11 = 1;
   # Check columns 13 to 42 (other diagnoses)
   for (i=13; i<42; i++) {
       if ($i == "G20") has_G20 = 1;
       if ($i == "E11") has_E11 = 1;
   }
   print $2, $3, has_G20, has_E11
}' G20_E11_NA.txt > updated_G20_E11.txt
```

Next, it's important that the data does not include any controls. The check step below can be used to ensure every entry in the ```updated_G20_E11.txt``` output is either a T2D or PD case or both.
```bash
# Ensure no entry has "0" and "0" (we only want diseases cases in the data)
awk '$3 == 0 && $4 == 0 {print $0}' updated_G20_E11.txt
```
But, this output has multiple entries of the same individuals with same or different diagnoses. This must be harmonized for easier handling downstream.
```bash
# Using the script below, filter out instances of duplicate individuals with multiple diagnoses for same condition
awk 'BEGIN {FS=OFS="\t"} 
NR==1 {print "ID", "PD", "T2D"; next}  
{
    id = $1  # IDs are now in column 1
    if (id in seen) {
        if ($3 == 1) has_G20[id] = 1
        if ($4 == 1) has_E11[id] = 1
    } else {
        seen[id] = 1
        has_G20[id] = $3
        has_E11[id] = $4
    }
} 
END {
    for (id in seen) {
        print id, has_G20[id], has_E11[id]
    }
}' updated_G20_E11.txt > individuals_with_PD_T2D.txt
```
Use the ID match document to update the individual IDs (to match the genotyping data)

```bash
# Update IDs
awk 'NR==FNR {a[$2]=$1; next} {if ($1 in a) $1=a[$1]; print}' ../../../ID_match_unix.txt individuals_with_PD_T2D.txt > PD_T2D_harmonized.txt

# Update covariates: merge the data in both the old covariates file and PD_T2D_harmonized.txt into one 
awk 'NR==FNR {t2d[$1] = $3; next} FNR==1 {print $0, "T2D"; next} {print $0, ($1 in t2d) ? t2d[$1] : "0"}' PD_T2D_harmonized.txt ../02_covariate_file/covariates.txt > new_covariates.txt
```

- This final output includes 948 individuals (832 PD, 55 T2D, and 61 both PD & T2D)
- At this point, it is necessary to compare with the questionnaire data to see if there are differences
- Questionnaire data has 935 cases altogether (867 PD, 54 T2D, and 62 both PD & T2D). 
- The national diagnostic data will be use as the updated version, specifically for studying the impcat of T2D on PD in the MPBC cohort.

```bash
# First, do some quality checks
# 1. check if all the IDs are in each file or not (not imprtant)
awk 'NR==FNR {ids[$1]; next} !($1 in ids) {print $1 " is missing from file B"}' PD_T2D_harmonized.txt covariates_with_T2D.txt | wc -l 
This should return only a list of healthy controls that should be in the dataset (916)
# 2. check how many PD cases were actually obtained from the new diagnostic data
awk '$2==1 {print $1, $2}' PD_T2D_harmonized.txt | grep '^GB' | wc -l
The number here is 893 (36 less than our original data)
# 3. Check how many T2D cases are there as well
awk '$3==1 {print $1, $3}' PD_T2D_harmonized.txt | grep '^GB' | wc -l
116 (32 less than we have in the original data)
 - Since the focus now is on T2D only, we need to update only info for T2D in the covariates file. 

# 4. check if the data for cases of PD and T2D match in the two datasets (not important)
awk 'NR==FNR {pd_status[$1]=$6; next} ($1 in pd_status) {print $1, pd_status[$1], $2}' covariates_with_T2D.txt PD_T2D_harmonized.txt | wc -l


# count the PD, T2D, PD + T2D, & healthy controls
```bash
[elugoh@sens12 08_T2D_PD]$ awk '$6==1 {print $20}' new_covariates.txt | sort | uniq -c
    867 0
     62 1
[elugoh@sens12 08_T2D_PD]$ awk '$6==0 {print $20}' new_covariates.txt | sort | uniq -c
    881 0
     54 1
```
# copy all four summary stats to the working dir
>sftp put -r Scores_PD-DM/ /home/elugoh/Documents/DataFolder/Elijah/processed_analysis_files/08_T2D_PD/
# update the formatting of the each summary stat file to tab-deliminated
sed 's/ /\t/g' Score_Ge.txt > tab_delim_Ge.txt
sed -i 's/ /\t/g' Score_Khera.txt 
sed -i 's/ /\t/g' Score_Lin.txt 
sed -i 's/ /\t/g' Score_Mars.txt
# next, to check for snp id and allele consistency and update the snp ids, the following python script is run
```flip_ref_alt.py``` using the SLURM job script, ```update_scores.sh```, on CSENS

- This step is repeated for all score files

4 summary stats are obtained from established studies and uploaded to the cluster. These are:
1. [Khera 2018](https://www.nature.com/articles/s41588-018-0183-z)
2. [Mars 2020](https://www.nature.com/articles/s41591-020-0800-0)
3. [Ge 2022](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01074-2), and 
4. [Lin 2023](https://www.sciencedirect.com/science/article/pii/S0048969723028747?via%3Dihub)


# Now, the PRS Calculations
```bash
mkdir PRS_T2D && cd PRS_T2D
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Khera_scores2.txt --out Khera_PRS
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Ge_scores2.txt --out Ge_PRS
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Mars_scores2.txt --out Mars_PRS
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Lin_scores2.txt --out Lin_PRS

# cp the PRS calculations and the new covariate file to the local computer and load them in R for further analysis and visualization

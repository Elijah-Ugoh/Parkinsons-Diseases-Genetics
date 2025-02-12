# Investigating the Impact of Type II Diabetes Mellitus (T2D) on Parkinson's Disease (PD) Risk in the MPBC Cohort
Steps:
- Data QC (extensive QC has been done on the MPBC cohort. So, no need for QC here, since the same dataset will be used in this pipeline).
- Create a new covariate file with additional information on T2D for all participants.
- To create the covariate file, we must first extract diagnostic information from government data to ensure we include only T2D (and not T1D) in the analysis.  
- Get summary stats from publicly availble GWAS data and update them to plink-ready format.
- Match the SNPs in the summary stat with those in our genotype data to ensure the alleles and ID consistency. 
- Run the PRS calculations using plink2. 
- Analyze the results using logistic regressisons in R. 

The steps below create an updated covariate file:
```bash

```

4 summary stats are obtained from established studies and uploaded to the cluster. These are:
1. [Khera 2018](https://www.nature.com/articles/s41588-018-0183-z)
2. [Mars 2020](https://www.nature.com/articles/s41591-020-0800-0)
3. [Ge 2022](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01074-2), and 
4. [Lin 2023](https://www.sciencedirect.com/science/article/pii/S0048969723028747?via%3Dihub)


# view the headers to see which columns must be extracted
less -S ../QuestionnaireData_N1864_FINAL_CLEANED_210621_nomed.csv | head -1 | tr "," "\n" | nl
# extract the diabetes information and save to a new file
cut -d ',' -f 88 ../QuestionnaireData_N1864_FINAL_CLEANED_210621_nomed.csv | sed 's/,/\t/g' | sed 's/"//g' > diabetes.txt
# update the IDs 
awk 'NR==FNR {a[$2]=$1; next} {if ($1 in a) $1=a[$1]; print}' ../ID_match_unix.txt diabetes.txt > diabetes_with_ID.txt 
# create an updated covariate file by merging in info from the diabetes file
header1=$(head -n 1 02_covariate_file/covariates.txt) # extract the headers first
header2=$(head -n 1 diabetes_with_ID.txt | awk '{print $NF}') 
echo "$header1 $header2" > covariates_with_T2D.txt #merge the headers
awk 'NR==FNR {a[$1] = $2; next} {if ($1 in a) print $0 "\t" a[$1]}' diabetes_with_ID.txt 02_covariate_file/covariates.txt >> covariates_with_T2D.txt #save the data to the same file as the headers
sed -i 's/ /\t/g' covariates_with_T2D.txt # ensure the dataset is properly tab-deliminated
# now move the all the files to the new working dir
mv covariates_with_T2D.txt diabetes.txt diabetes_with_ID.txt 08_T2D_PD/
# cd into it
cd 08_T2D_PD/
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

Next, since the questionnaire data does not contain info as to whether the diabetes diagnosis is Type 1 or Type 11, we use supplementary data from General diagnositic information of the patients
# the individuals with either or both E11 (T2D) of G20 (PD) are filtered out to a separate file
(head -n 1 UT_R_PAR_OV_14691_2021.txt && grep -E 'E11|G20' UT_R_PAR_OV_14691_2021.txt) > G20_E11.txt
# assign NAs to first columns with missing entries
awk 'BEGIN {FS=OFS="\t"} {if ($1 == "") $1="NA"; print}' G20_E11.txt > G20_E11_NA.txt
# print the relevant columns and redirect the output to a different file
```bash
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
# Ensure no entry has "0" and "0" (no controls, only cases)
awk '$3 == 0 && $4 == 0 {print $0}' updated_G20_E11.txt

# Next, filter out the duplicates with multiple diagnoses for same condition
bash filter_dupliates.sh  
# Update the IDs
awk 'NR==FNR {a[$2]=$1; next} {if ($1 in a) $1=a[$1]; print}' ../ID_match_unix.txt individuals_with_PD_T2D.txt > PD_T2D_harmonized.txt
This final file includes 948 individuals with either PD, T2D or both in the MPBC cohort
We need to compare with the covariates data to ensure they are the same, because we initially only have 935 (928 with PD only & 77 with T2D only) in that
The national diagnostic data will be use as the updated version

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

# Now, merge the data from both file into one (this is the next step after creating the (PD_T2D_harmonized.txt) file. 
awk 'NR==FNR {t2d[$1] = $3; next} FNR==1 {print $0, "T2D"; next} {print $0, ($1 in t2d) ? t2d[$1] : "0"}' PD_T2D_harmonized.txt covariates.txt > new_covariates.txt

# count the PD, T2D, PD + T2D, & healthy controls
```bash
[elugoh@sens12 08_T2D_PD]$ awk '$6==1 {print $20}' new_covariates.txt | sort | uniq -c
    867 0
     62 1
[elugoh@sens12 08_T2D_PD]$ awk '$6==0 {print $20}' new_covariates.txt | sort | uniq -c
    881 0
     54 1
```
# Now, the PRS Calculations
```bash
mkdir PRS_T2D && cd PRS_T2D
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Khera_scores2.txt --out Khera_PRS
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Ge_scores2.txt --out Ge_PRS
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Mars_scores2.txt --out Mars_PRS
plink --bfile ../../../MPBC_HRC_Rsq03_updated --score ../Scores_PD-DM/updated_Lin_scores2.txt --out Lin_PRS

# cp the PRS calculations and the new covariate file to the local computer and load them in R for further analysis and visualization

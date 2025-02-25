# Stratification of Polygenic Risk Scores (PRS) by Lifestyle Factors in Parkinson’s Disease
## Overview

Here the PRS scores obtianed in ```5_PRS_Computation.md``` is stratified by different lifestyle/environmental factors, whi h in this case are smoking, caffeine consumption, and pesticide exposure, to examine their interaction with Parkinson’s disease (PD). The goal is to assess whether PRS has a synergistic effect on PD with these factors, as well as assess how PRS distributions vary across these factors.

Dataset Description
The dataset contains the following relevant columns:
- PD is categorized as a binary variable (1 = Parkinson’s case, 0 = control).
- Smoking_Ever: Whether the individual has ever smoked (YES/NO).
- Snus_Ever: Whether the individual has used Snus (tobacco pouch) (YES/NO).
- Caffeine_Ever: Whether the individual has ever drank coffee, black tea, and green tea at any time below the age of 40, between 41-64 years, and above 64 years. (YES/NO).
- Caffeine_Consumption_Level: Categorized as LOW, MEDIUM, or HIGH.
- Pesticide_Exposure: Whether the indiviadual has been exposed to/used pesticides before or not (YES/NO).
- Covariates for logistic regressions include age, sex, education level, and first 5 principal components.

Stratification Steps
The data is prepared as follows: 
- Exclude individuals with missing data for any lifestyle factor and extract each lifestyle factor data in separate files.
- Compute mean and standard deviation of PRS for each stratified group.
- Run logistic regression (adjusted for covariates) to assess additive effects of PRS and each factor.
- Visualize PRS distributions using violin plots and interaction plots.
- Test for statistical interactions between PRS and lifestyle factors.

## Extraction of Lifestyle Factor Data from Questionnaire
```bash
# standing in the main wd, create a new directory
mkdir processed_analysis_files/07_stratification
cut -d ',' -f 2,4,6,10,24-51,88,179,181 QuestionnaireData_N1864_FINAL_CLEANED_210621_nomed.csv | sed 's/,/\t/g' > processed_analysis_files/07_stratification/env_factors.txt
# Update the individual IDs
awk 'NR==FNR {a[$2]=$1; next} {if ($1 in a) $1=a[$1]; print}' /ID_match_unix.txt processed_analysis_files/07_stratification/env_factors.txt > processed_analysis_files/07_stratification/new_env_factors.txt
```
PRS Scores obtained from the 90 variants identified by Nalls et al is used. These scores are merged into the data just extracted above.  
```bash
# now, merge the scores column from the score file with the rest of the environment data just created above
```bash
header1=$(head -n 1 new_env_factors.txt) # extract the headers first
header=$(../04_PRS/PRS_PD_MPBC.profile | awk '{print $NF}')
echo "$header1 $header2" > merged_score_and_env.txt #merge and save to a new file

# add the data to the headers (the process here involves using the FIDs to match the scores to the right individuals)
awk 'NR==FNR {a[$1] = $NF; next} {if ($1 in a) print $0 "\t" a[$1]}' ../04_PRS/PRS_PD_MPBC.profile new_env_factors.txt >> merged_score_and_env.txt
sed -i 's/ /\t/g' merged_score_and_env.txt # re-format and make data tab-deliminated
```
Now, extract the three environmental lifestyle categories and create different files
```bash
cut -f 1-4,9,34-36 merged_score_and_env.txt > smoking.txt
cut -f 1-4,15,34-36 merged_score_and_env.txt > snus.txt
cut -f 1-5,34-36 merged_score_and_env.txt > pesticides.txt
cut -f 1-4,24-32,34-36 merged_score_and_env.txt > caffeine.txt 
```
Add a new column for either binary or categorical traits
NB: The following blocks of code also take into account missing entries "NA" and the delimiter, "\t", in the files. 

- Smoking
```bash
awk -F"\t" '
BEGIN {OFS="\t"}
NR==1 {print $0, "Smoking_Ever"}
NR>1 {
    if ($5 !="NA" && $5 > 0)
        print $0, "YES";
    else if ($5 !="NA")
        print $0, "NO";
}' smoking.txt > smoking_ever_never.txt
```
- Using snus
```bash
awk -F"\t" '
BEGIN {OFS="\t"}
NR==1 {print $0, "Snus_Ever"}
NR>1 {
    if ($5 !="NA" && $5 > 0)
        print $0, "YES";
    else if ($5 !="NA")
        print $0, "NO";
}' snus.txt > snus_ever_never.txt
```
- Pesticide exposure
```bash
awk -F"\t" '
BEGIN {OFS="\t"}
NR==1 {print $0, "Pesticides_Ever"}
NR>1 {
    if ($5 !="NA" && $5 ==2)
        print $0, "YES";
    else if ($5 !="NA" $$ $5 ==1)
        print $0, "NO";
}' pesticides.txt > pesticides_ever_never.txt
```
- Caffeine consumption (no cup vs at leats 1 cup/day)
```bash
awk -F"\t" '
BEGIN {OFS="\t"}
NR==1 {print $0, "Caffeine_Ever"}
NR>1 {
    if ($5 !="NA" && $5 > 1 || $6 !="NA" && $6 > 1 || $7 !="NA" && $7 > 1 || $8 !="NA" && $8 > 1 || $9 !="NA" && $9 > 1 || $10 !="NA" && $10 > 1 || $11 !="NA" && $11 > 1 || $12 !="NA" && $12 > 1 || $13 !="NA" && $13 > 1)
        print $0, "YES";
    else if ($5 !="NA" $$ $5 ==1)
        print $0, "NO";
}' caffeine.txt > caffeine_ever_never.txt
```
- Caffeine consumption - low (1-2 cups/day) vs high (3-5 cups/day)
```bash
awk -F"\t" '
BEGIN {OFS="\t"}
NR==1 {print $0, "Caffeine_Consumption_Level"}
NR>1 {
    if ($5 !="NA" && $5 > 2 || $6 !="NA" && $6 > 2 || $7 !="NA" && $7 > 2 || $8 !="NA" && $8 > 2 || $9 !="NA" && $9 > 2 || $10 !="NA" && $10 > 2 || $11 !="NA" && $11 > 2 || $12 != "NA" && $12 > 2 || $13 !="NA" && $13 > 2)
        print $0, "HIGH";
    else if ($5 !="NA" && $5 < 3 || $6 !="NA" && $6 < 3 || $7 !="NA" && $7 < 3 || $8 !="NA" && $8 < 3 || $9 !="NA" && $9 < 3 || $10 !="NA" && $10 < 3 || $11 !="NA" && $11 < 3 || $12 !="NA" && $12 < 3 || $13 !="NA" && $13 < 3)
        print $0, "LOW";
}' caffeine.txt > caffeine_low_vs_high.txt
```
- Next, copy the generated files to local and read them into R and run the 5 logistic regresions written in the ```stratification_analysis.R``` script.

NB: All statistical analysis and visualisations for the stratification is done using the ```stratification_analysis.R``` script. 

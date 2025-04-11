# Investigating the Associations Between Polygenic Risk Scores, Environmental Exposure, and Type II Diabetes in Parkinson’s Disease

## Overview
This repository contains all the scripts, resources, and analysis workflow used in the study titled:  
*"Investigating the Associations Between Polygenic Risk Scores, Environmental Exposure, and Type II Diabetes in Parkinson’s Disease"*.

The study aimed to evaluate gene-environment interaction in Parkinson's disease (PD) risk using polygenic risk scores (PRS) and environmental and lifestyle exposure data (smoking, caffeine, snus, pesticide), and the potential genetic contribution of Type II Diabetes (T2D) to PD risk.



## Getting Started
The complete analysis for this study can be followed in this order:

1. [Setup and Installations](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/1_Setup_and_Installations.md)
2. [SNPs Extraction and Formatting](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/2_extract_90_risk_variants.md)
3. [Extraction of Covariate Data from Questionnaire](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/3_creating_covariates.md)
4. [Genotyping Data Quality Control](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/4_quality_control.md)
5. [PRS Computations](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/5_PRS_Computation.md)
6. [Stratification Analyses](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/6_stratification.md)
7. [Investigating T2D Impact on PD](https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics/blob/main/README.md).

NB: Some code for analysis were inspired by the pipeline in the [Global Parkinson's Disease Genetics programme](https://github.com/GP2-TNC-WG/GP2-Beginner-Bioinformatics-for-PD-Genetics/blob/master/README.md).

## To Reproduce this Analysis
1. Clone the repository:
```bash
git clone https://github.com/Elijah-Ugoh/Parkinsons-Diseases-Genetics.git
```
## Acknowledgements
This project was completed at and with help from the Translational Neurogenetics Unit, Faculty of Medicine, Lund University, Sweden.

All analyses utilized questionnaire and genotyping data from 1902 participants in the [MultiPark Biobank Sample Collection (MPBC) ](https://www.multipark.lu.se/infrastructures/biobank-platform) study cohort, including cases and controls. T2D analysis utilized data from [Socialstyrelsen](https://www.socialstyrelsen.se/en/statistics-and-data/registers/national-patient-register/).

### Citation
If you use any part of this code or workflow, please cite the repository and the related manuscript (when published).

### Contact
For questions or contributions, please feel free to open an issue or submit a pull request on GitHub.

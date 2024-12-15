"""
#!/usr/bin/env python3
author: Elijah C. Ugoh
Date: 11/19/2024

Use this script to update the PHENO column in the .bim file befor running GWAS in Plink
"""

import sys

def load_fam_file(fam_file):
    """Load the .fam file and return a dictionary with FID as key and phenotype as value."""
    fam_data = {}
    with open(fam_file, "r") as infile:
        for line in infile:
            fields = line.strip().split()
            FID = fields[0]
            phenotype = fields[5]
            fam_data[FID] = phenotype
    return fam_data

def update_covariate_file(covariate_file, fam_data, output_file):
    #Update the covariate file with phenotypes from the fam data
    with open(covariate_file, "r") as covariates, open(output_file, "w") as outfile:
        # Process the header
        header = covariates.readline().strip()
        outfile.write(header + "\n")
        
        for line in covariates:
            fields = line.strip().split()
            FID = fields[0]
            
            if FID in fam_data:
                fields[-1] = fam_data[FID]  # Update the phenotype
            else:
                print(f"Warning: Individual {FID} in covariate file not found in fam file. Skipping.")
            
            outfile.write(" ".join(fields) + "\n")

# Execution block
if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python update_pheno.py <fam_file> <covariate_file> [output_file]")
        sys.exit(1)

    fam_file = sys.argv[1]
    covariate_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) == 4 else "updated_covariates.txt"

    # Load data and update the covariate file
    fam_data = load_fam_file(fam_file)
    update_covariate_file(covariate_file, fam_data, output_file)

    print(f"Updated covariate file saved as {output_file}")
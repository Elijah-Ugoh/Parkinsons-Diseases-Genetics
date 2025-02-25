"""
- This script only works if the bim file has the format chr:pos:ref:alt
- This script checks for exact SNP matches between a GWAS PRS score file and a .bim file.
- If a match is found (based on chr, pos, effect allele matching the ref allele in bim, and/or a flipped ref/alt alleles in the bim file), the SNP ID in the score file is updated by replacing it with the corresponding SNP ID from the .bim file.
- If no match is found, the SNP is written to a separate file for manual review.
- No SNPs are lost; all SNPs from the score file are accounted for.
"""

import time
import pandas as pd

# start time tracking
start_time = time.time()

print("Processing. Please grab a cup of coffee or water while you wait...")

# Load the PRS score file
score_file = "Score_Mars.txt"
score_df = pd.read_csv(score_file, sep="\t", dtype={"chr_name": str, "chr_position": int})

# Load bim file
bim_file = "MPBC_HRC_Rsq03_updated.bim"
bim_df = pd.read_csv(bim_file, sep="\t", header=None, 
                     names=["bim_chr", "bim_snp_id", "chr_dist", "pos", "bim_ref", "bim_alt"],
                     dtype={"bim_chr": str, "pos": int})

# Create a lookup dictionary for quick matching (key = (chr, pos), value = (snp_id, ref, alt))
bim_lookup = {(row["bim_chr"], row["pos"]): (row["bim_snp_id"], row["bim_ref"], row["bim_alt"])
              for _, row in bim_df.iterrows()}

# Prepare lists for output files
updated_snps = []
unmatched_snps = []

# Iterate over the score file
for _, row in score_df.iterrows():
    key = (row["chr_name"], row["chr_position"])
    
    if key in bim_lookup:
        bim_snp_id, bim_ref, bim_alt = bim_lookup[key]
        
        if row["effect_allele"] == bim_ref and row["other_allele"] == bim_alt:  # Exact match found
            updated_snps.append([bim_snp_id, row["effect_allele"], row["effect_weight"]])
        elif row["effect_allele"] == bim_alt and row["other_allele"] == bim_ref:  # swapped alleles match
            updated_snps.append([bim_snp_id, row["effect_allele"], row["effect_weight"]])
        else:  # Position matched but alleles didn't
            unmatched_snps.append([row["rsID"], row["chr_name"], row["chr_position"], 
                                   row["effect_allele"], row["other_allele"], row["effect_weight"]])
    else:
        unmatched_snps.append([row["rsID"], row["chr_name"], row["chr_position"], 
                               row["effect_allele"], row["other_allele"], row["effect_weight"]])

# Convert lists to DataFrames and save to files
updated_score_df = pd.DataFrame(updated_snps, columns=["rsID", "effect_allele", "effect_weight"])
unmatched_snps_df = pd.DataFrame(unmatched_snps, columns=["rsID", "chr_name", "chr_position", 
                                                           "effect_allele", "other_allele", "effect_weight"])

# Save the corrected PRS score file **without headers**
updated_score_df.to_csv("updated_Mars_scores2.txt", sep="\t", index=False, header=False)

# Save unmatched SNPs **with headers** for manual review
unmatched_snps_df.to_csv("unmatched_Mars_snps2.txt", sep="\t", index=False, header=True)

# End time tracking
end_time = time.time()
# calculate and print execution time
execution_time = end_time - start_time
print(f"Script Run time:  {execution_time:.2f} seconds")
print(f"Updated PRS score file saved as: updated_Mars_scores2.txt\n"
      f"Unmatched SNPs saved to: unmatched_Mars_snps2.txt")

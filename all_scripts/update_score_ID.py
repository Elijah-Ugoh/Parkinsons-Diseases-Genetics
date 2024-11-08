"""
#!/usr/bin/env python3

This script updates the SNPIDs in your beta score file to match that of the .bim file where the .bim file has alt and ref allele attached to the SNP IDs.

How It Works:
- load_bim_file: Reads the .bim file and creates a dictionary where the key is Chromosome:Position and the value is the full SNP identifier.
- update_score_file: Reads the score file, matches each Chromosome:Position with the .bim dictionary, updates the SNP IDs, and writes the result to a new file.
- Handles exceptions and outputs informative error messages.

Usage: 
./update_score.py test.bim META5_GRS_chr_bp.txt updated_score_file.txt

"""

import sys 

# Function to load and process the .bim file and extract the SNP IDs as chr:pos only
def load_bim_file(bim_file):
    bim_data = {} 
    with open(bim_file, 'r') as bim:
        for line in bim:  # Loop through each line in the .bim file
            cols = line.strip().split()  
            # Combine the chromosome and position as the key
            snp_id = ':'.join(cols[1].split(':')[:2])  
            bim_data[snp_id] = cols[1]  
    return bim_data  # Return the dictionary containing SNP IDs

# Function to update the score file by matching IDs from the .bim file
def update_score_file(score_file, bim_data, output_file):
    with open(score_file, 'r') as score, open(output_file, 'w') as out:
        # Open the score file for reading and output file for writing
        for line in score:  
            cols = line.strip().split()  
            snp_id = cols[0] 
            # If the SNP ID exists in the .bim dictionary, replace it with the full ID
            if snp_id in bim_data:
                cols[0] = bim_data[snp_id] 
            else:
                # If the SNP ID is not found in the .bim file, print a warning
                print(f"Warning: SNP ID {snp_id} not found in .bim file, keeping original ID")
            out.write("\t".join(cols) + "\n")  # Write the updated line to the output file

# Execution block
if __name__ == "__main__":
    # Ensure the correct number of command-line arguments are provided
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python script.py <bim_file> <score_file> [output_file]")
        sys.exit(1)  # Exit if the arguments are incorrect

    bim_file = sys.argv[1]  
    score_file = sys.argv[2]
    # If an output file is specified, use it; otherwise, overwrite the score file
    output_file = sys.argv[3] if len(sys.argv) == 4 else score_file

    bim_data = load_bim_file(bim_file)  
    update_score_file(score_file, bim_data, output_file)

    print(f"Updated score file saved as {output_file}")  # Final message to notify the user of the output

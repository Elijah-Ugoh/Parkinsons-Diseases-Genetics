import pandas as pd

# Load the reference file
reference_df = pd.read_csv('updated_score_file.txt, sep='\t')
# Adjust column names as needed; here we assume 'ID' is rsID, 'CHR' is chromosome, 'BP' is position, 'Allele1' and 'Allele2' are alleles
reference_df['chr_pos_ref_alt'] = (
    reference_df['CHR'].astype(str) + ":" +
    reference_df['BP'].astype(str) + ":" +
    reference_df['Allele1'].astype(str).str.upper() + ":" +
    reference_df['Allele2'].astype(str).str.upper()
)
# Create a dictionary for fast lookup of rsID to chr:pos:ref:alt
rsid_to_chrpos = dict(zip(reference_df['ID'], reference_df['chr_pos_ref_alt']))

# Load the score file
score_df = pd.read_csv('new_rsIDs.txt', sep='\t', header=None, names=['SNP', 'A2', 'b'])

# Replace rsID with chr:pos:ref:alt format using the dictionary
score_df['chr_pos_ref_alt'] = score_df['SNP'].map(rsid_to_chrpos)

#drop rows with misisng  SNP IDs
score_df = score_df.dropna(subset=['chr_pos_ref_alt'])

# Drop the old rsID column if no longer needed and save
score_df.drop(columns=['SNP'], inplace=True)

#Reorder columns so chr_pos_ref_alt comes first, then EffectAllele, and then Beta
score_df = score_df[['chr_pos_ref_alt','A2','b']]

#save to file 
score_df.to_csv('cleaned_score_file.txt', sep= '\t', index=False, header=False)

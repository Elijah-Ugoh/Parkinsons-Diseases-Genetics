# Creating a Score File for PRS Calculation
The summmary statistics for only 90 significant PD risk loci (SNPs) are obtained from the [Nalls study supplementary appendix](https://pmc.ncbi.nlm.nih.gov/articles/instance/8422160/bin/NIHMS1735188-supplement-mmc1.pdf). Click on the [additional supplement](https://drive.google.com/file/d/1VUU_-tYl-ew08vupEVRuNEPDpTpWWom_/view?usp=sharing) link to download the file.

- The file of interest is ```Table S2. Detailed summary statistics on all nominated risk variants, known and novel_.xlsx```. Download to a preferred local dir.
- Next, we filter out the columns of interest, including the SNP IDs, alleles, and beta score columns from the  meta analysis.

```bash
mkdir cosmos-home  # create local working dir and save the downloaded file here
# install command line tool for converting xlsx files to csv files on Linux
sudo apt install gnumeric
# Convert the spreadsheet to CSV
ssconvert Table\ S2.\ Detailed\ summary\ statistics\ on\ all\ nominated\ risk\ variants\,\ known\ and\ novel_.xlsx  meta_analyis_90.csv

less meta_analyis_90.csv  # View the CSV content

# extract the columns of interest and save to a new file
awk -F, '{print $1,$2,$3,$6,$7,$9,$10}' meta_analyis_90.csv > meta_analyis_90.txt

# Re-format the SNP IDs to align with the .bim file in the genotype data (chr:pos:ref:alt)
awk 'NR>1 {print $2":"$3":"toupper($4)":"toupper($5)"\t"toupper($4)"\t"toupper($5)"\t"$7}' meta_analyis_90.txt > meta_90.txt

# remove all columns and rows containing information we don't want. Only first 90 lines needed. 
grep -E '^[0-9]+:[0-9]+' meta_90.txt | head -n 90 > metaanalysis90.txt

# Command explanation:
"^[0-9]+:[0-9]+" matches the chromosome format
"head -n 90" prints first 90 lines.

# finally, check the content of the cleaned file and count the number of rows
less metaanalyis90.txt
wc -l metaanalyis90.txt
```

Move the cleaned score file from local working dir to working dir on the cluster
```bash
sftp elugoh@cs-diode.lunarc.lu.se
put cosmos-home/metaanalyis90.txt /home/elugoh/Documents/DataFolder/Elijah
```


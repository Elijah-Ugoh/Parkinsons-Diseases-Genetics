This analysis is first run using the beta score for only 90 risk SNPs associated with PD and obtained from: [Nalls et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/instance/8422160/bin/NIHMS1735188-supplement-mmc1.pdf)

- First, we filter out the columns of interests along with the meta analysis beta score of the identified risk variants. 
- The file of interest is ```Table S2. Detailed summary statistics on all nominated risk variants, known and novel_.xlsx```.

```bash 
# install comand line tool for converting xlsx files to csv files on Linux
sudo apt install gnumeric
# Convert the spreadsheet to CSV
ssconvert ../Table\ S2.\ Detailed\ summary\ statistics\ on\ all\ nominated\ risk\ variants\,\ known\ and\ novel_.xlsx  meta_analyis_90.csv

less meta_analyis_90.csv  # View the CSV content

# extract the columns of interest and save to a new file
awk -F, 'NR>1 {print $1":"$2"\t"toupper($3)"\t"$4}' meta_analyis_90.csv > meta_analyis_90.txt

# remove all columns and rows containing information we don't want
grep -E '^[0-9]+:[0-9]+\s+[ACGT]\s+[-+]?[0-9]+([.][0-9]+)?' meta_analyis_90.txt > metaanalyis90.txt

# Command explanation:
^[0-9]+:[0-9]+ matches the chromosome format.
\s+[ACGT]\s+ matches a single nucleotide (A, C, G, or T).
[-+]?[0-9]+([.][0-9]+)? matches the beta values.

# finally, check the content of the cleaned file and count the number of rows
less metaanalyis90.txt
wc -l metaanalyis90.txt
```

Move the beta score file from local working dir to remote working dir on the cluster
```bash
sftp elugoh@cs-diode.lunarc.lu.se
put cosmos-home/metaanalyis90.txt /home/elugoh/Documents/DataFolder/Elijah
```


# Quality Control and Imputation
To avoid spurious associations and bias. QC is performed at both sample and variant level on the genotyping data.

The file needed are the raw Plink binary files (.fam, .bim, and .bed files) for pd patienst and controls. PLINK and GCTA are needed here.

```bash
mkdir processed_analysis_file/01_QC

# check if plink is installed and load it
module spider plink
module spider PLINK/2.00a3.7
module load GCC/11.3.0  OpenMPI/4.1.4 # load dependencies
module load PLINK/2.00a3.7 # then load plink
```
## Genotyping Daat QC
### 1. Filter Based on Genotyping Call Rate per Sample
```bash
plink --bfile MPBC --missing --out processed_analysis_files/01_QC/call_rates
```

#### Remove Call Rate Outliers
- Based on the results from the ```processed_analysis_files/1_QC/call_rates.imiss``` file, samples that are not genotyped at the 95% of the tested variance are removed using the ```--mind 0.05``` flag.
- No individual was removed as all passed the filtering with good genotyping call rate of 0.97.

```bash
plink --bfile MPBC --mind 0.05 --make-bed --out processed_analysis_files/01_QC/MPBC_call_rates
```

NB: No individual was removed due to missing genotype data 

### 2. Next, Prune the Data for Heterozygosity and Linkage Disequilibrium
```bash
mkdir processed_analysis_file/01_QC/prunning
plink --bfile MPBC --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out processed_analysis_files/01_QC/prunning/prunning
plink --bfile MPBC --extract processed_analysis_files/05_QC/prunning/prunning.prune.in --make-bed --out processed_analysis_files/01_QC/prunning/prunned_data
cd processed_analysis_files/01_QC/prunning/ && plink --bfile prunned_data --het --out prunnedHet

# save the outliers
awk '{if ($6 <= -0.15) print $0}' prunnedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0}' prunnedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > Heterozygosity_Outliers.txt
cut -f 1,2 Heterozygosity_Outliers.txt > all_outliers.txt

# remove the outlying individuals
plink --bfile ../../01_QC/MPBC_call_rates --remove all_outliers.txt --make-bed --out after_het_call_rate
```

NB: Only 1 individual was removed after filtering for heterozygosity. 

### 3. Proceed to Sample-level and Sex QC 
```bash
plink --bfile after_het_call_rate --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
```

NB: No sex check problem found. So, no idividuals was filtered out at this stage. 

```bash
# Make binary files to proceed
plink --bfile after_het_call_rate --make-bed --out after_gender
```

### 4. Do Ancestry QC
- The samples in the plink binary files are merged with a reference panel, i.e., [HAPMAP reference populations](http://www.hapmap.org/). This way, all samples in the study that are not from the ancestry of interest (Europe) can be filtered out. 
- This comparison with HAPMAP is based on the number of SNPs in the input dataset that overlap with the HAPMAP data. 

```bash
# move up one directory and make a new directory
cd .. && mkdir ancestry

# Attempt to merge the target dataset with the HapMap dataset to identify common SNPs.
plink --bfile prunning/after_gender --bmerge ../../HAPMAP_files/HAPMAP_hg19_new --out ancestry/hapmap3_bin_snplis --make-bed

# Correct mismatched alleles (flipped SNPs) from the failed merge using the missnp file.
plink --bfile prunning/after_gender --flip ancestry/hapmap3_bin_snplis-merge.missnp --make-bed --out ancestry/after_gender3 

# Retry merging the corrected dataset with the HapMap dataset after flipping SNPs.
plink --bfile ancestry/after_gender3 --bmerge ../../HAPMAP_files/HAPMAP_hg19_new --out ancestry/hapmap3_bin_snplis --make-bed 

# Remove SNPs that still fail to merge after flipping, creating a new dataset.
plink --bfile ancestry/after_gender3 --exclude ancestry/hapmap3_bin_snplis-merge.missnp --out ancestry/after_gender4 --make-bed

# Perform the final merge after excluding problematic SNPs.
plink --bfile ancestry/after_gender4 --bmerge ../../HAPMAP_files/HAPMAP_hg19_new --out ancestry/hapmap3_bin_snplis --make-bed

# Perform SNP pruning by excluding SNPs with high missingness, low MAF, and high LD.
plink --bfile ancestry/hapmap3_bin_snplis --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out ancestry/pruning

# Create a new dataset including only the pruned SNPs that passed the filtering criteria.
plink --bfile ancestry/hapmap3_bin_snplis --extract ancestry/pruning.prune.in --make-bed --out ancestry/pruned_data

# Calculate heterozygosity rates to identify potential sample quality issues.
plink --bfile ancestry/pruned_data --het --out ancestry/prunedHet

# Perform principal component analysis (PCA) using pruned SNPs to capture population structure.
plink --bfile ancestry/pruned_data --geno 0.01 --out ancestry/pca --make-bed --pca 4
```

NB: Duplicate SNPs with multiple positions in the hapmap data and target data were removed. No indiviadual was removed. 

#### Calculate PCA
```bash
cd ancestry 
mkdir pca_filtering
grep "EUROPE" pca.eigenvec > pca_filtering/eur.txt
grep "ASIA" pca.eigenvec > pca_filtering/asia.txt
grep "AFRICA" pca.eigenvec > pca_filtering/afri.txt
grep -v -f pca_filtering/eur.txt pca.eigenvec | grep -v -f pca_filtering/asia.txt | grep -v -f pca_filtering/afri.txt > pca_filtering/new_samples.txt
cut -d " " -f 3 ../prunning/after_gender.fam > pca_filtering/new_samples_add.txt
```

- Add unique identifiers to the ID's in each sample group, then merge them again to perform the PC
- "0" for the MPBC and European population, "1" for the european population, "2" and "3", for the Asian and African population, respectively
- This is necessary for separating the different data groups for PCA calculation.

```bash
cd pca_filtering
paste new_samples_add.txt new_samples.txt > new_samples2.txt
awk '{print 1, $0}' eur.txt > euro.txt
awk '{print 2, $0}' asia.txt > asiao.txt
awk '{print 3, $0}' afri.txt > afrio.txt 

# concatenate all the files into one
cat new_samples2.txt euro.txt asiao.txt afrio.txt > pca.eigenvec2

# plot the ancestry clusters
module load GCC/12.3.0
module load R/4.3.2
R < ../../../../HAPMAP_files/PCA_in_R.R --no-save
```
 
Finally, on acestry filtering, we keep the samples that cluster with Europeans in the HAPMAP data and make genotype binary files from those.

```bash
plink --bfile ../../prunning/after_gender --keep PCA_filtered_europeans.txt --make-bed --out after_gender_heterozyg_hapmap

# save all the outliers to a different file
cat PCA_filtered_asians.txt PCA_filtered_africans.txt PCA_filtered_mixed.txt > hapmap_outliers.txt 
```

After the PCA calculation, 11 individuals that didn't cluster well with the European ancestry were removed. 1890 individual are left. 

### 5. Filter based on Relatedness
```bash
cd ../../ && mkdir relatedness
plink --bfile ancestry/pca_filtering/after_gender_heterozyg_hapmap --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out relatedness/pruning
plink --bfile ancestry/pca_filtering/after_gender_heterozyg_hapmap --extract relatedness/pruning.prune.in --make-bed --out relatedness/pruned_data
plink --bfile relatedness/pruned_data --het --out relatedness/prunedHet
```
Run GCTA to Calculate Genetic Relatedness Matrix
- The Linux ```gcta-1.94.1-linux-kernel-3-x86_64.zip``` version is used and downlaoded from the [GCTA website](https://yanglab.westlake.edu.cn/software/gcta/#Download).

```bash
cd relatedness
# run gcta to compute genetic relatedness matrix using 5 threads
./gcta64 --bfile pruned_data --make-grm --out GRM --autosome --maf 0.05 --thread-num 5
# run the relatedness threshold cutoff
./gcta64 --grm-cutoff 0.125 --grm GRM --out GRM_0125 --make-grm --thread-num 5
```
NB: After pruning the GRM, 22 individuals were removed, leaving 1868 left in the dataset. 

Now, use the GCTA output files to complete the relatedness computation in plink
```bash
plink --bfile ../ancestry/pca_filtering/after_gender_heterozyg_hapmap --keep GRM_0125.grm.id --make-bed --out after_gender_heterozyg_hapmap_pihat
```

- The related individuals can also be filtered out directly using IBD (Identity By Descent) in plink. But, this filters out alternate individuals among related pairs. 
  
```bash
mkdir using_plink && cd using_plink

# perform linkage disequilibrium (LD) pruning to create a list of independent SNPs
plink --bfile ../../ancestry/pca_filtering/after_gender_heterozyg_hapmap --indep-pairwise 50 5 0.05 --out PPMI.pruned
# use only the independent SNPs listed in PPMI.pruned.prune.in for the IBD calculation
plink --bfile ../../ancestry/pca_filtering/after_gender_heterozyg_hapmap --extract PPMI.pruned.prune.in --genome --min 0.1 --out PPMI.IBD
# from the PPMI.IBD.genome file, which has information about related pairs, we extract the related pairs of individuals using an IBD threshold of 0.125 
awk '$10 > 0.125 {print $1, $2}' PPMI.IBD.genome | sort | uniq > IBD_remove.txt
# remove individuals listed in IBD_remove.txt
plink --bfile ../../ancestry/pca_filtering/after_gender_heterozyg_hapmap --remove IBD_remove.txt --make-bed --out after_gender_heterozyg_hapmap.pihat
# re-run IBD calculation on the new dataset to confirm that the remaining individuals are unrelated
plink --bfile after_gender_heterozyg_hapmap.pihat --extract PPMI.pruned.prune.in --genome --min 0.1 --out PPMI.IBD_check
```

NB: Here, 26 related individuals have been removed, as saved in the ```IBD_remove.txt```file above.

- Save IDs before and after the relatedness filteration
```bash
cut -f 1,2 ../../ancestry/pca_filtering/after_gender_heterozyg_hapmap.fam > IDs_before_relatedness_filter.txt
cut -f 1,2 after_gender_heterozyg_hapmap_pihat.fam > IDs_after_relatedness_filter.txt¢
```
- At the end of this filteration, plink removed 26 individual for relatedness, retaining only 1864 unrelated individuals. This is done for comparison only. The GCTA output is used subsequently.

### 6. Filter Based on Missingness per Variant
- This segment filters genetic variants based on missingness per variant, which helps to ensure data quality by removing SNPs with high levels of missing data. 

NB: Continue the remaining QC steps with outputs from GCTA. For reference, I have also provided the last two steps using the relatedness filtering output from plink above. 

```bash
# With GCTA Output
cd ..
mkdir missingness_per_variant && cd missingness_per_variant
# filter out SNPs with more than 5% missing genotype data (--geno 0.05) and create a new binary dataset 
plink --bfile ../relatedness/after_gender_heterozyg_hapmap.pihat --make-bed --out after_gender_heterozyg_pihat_mind --geno 0.05
# extract information on the SNPs removed by the --geno filter from the PLINK log file and save it. This document can be useful for tracking variants filtered due to high missingness
grep "(--geno)" after_gender_heterozyg_pihat_mind.log > MISSINGNESS_SNPS.txt
# test each SNP for differences in missing genotype rates between cases and controls, which can reveal genotyping quality issues specific to either group (Missingness by case control P > 1E-4)
plink --bfile after_gender_heterozyg_pihat_mind --test-missing --out missing_snps
awk '{if ($5 <= 0.0001) print $2}' missing_snps.missing > missing_snps_1E4.txt
# remove SNPs with significant case-control missingness difference
plink --bfile after_gender_heterozyg_pihat_mind --exclude missing_snps_1E4.txt --make-bed --out after_gender_heterozyg_pihat_mind_missing1
# create a unique, sorted list of excluded variants
sort -u missing_snps_1E4.txt > VARIANT_TEST_MISSING_SNPS.txt
```

```bash
# With Plink Output
cd ..
mkdir missingness_per_variant && cd missingness_per_variant
# filter out SNPs with more than 5% missing genotype data (--geno 0.05) and create a new binary dataset 
plink --bfile ../relatedness/using_plink/after_gender_heterozyg_hapmap.pihat --make-bed --out after_gender_heterozyg_pihat_mind --geno 0.05
# extract information on the SNPs removed by the --geno filter from the PLINK log file and save it This document can be useful for tracking variants filtered due to high missingness
grep "(--geno)" after_gender_heterozyg_pihat_mind.log > MISSINGNESS_SNPS.txt
# test each SNP for differences in missing genotype rates between cases and controls, which can reveal genotyping quality issues specific to either group (Missingness by case control P > 1E-4)
plink --bfile after_gender_heterozyg_pihat_mind --test-missing --out missing_snps
awk '{if ($5 <= 0.0001) print $2}' missing_snps.missing > missing_snps_1E4.txt
# remove SNPs with significant case-control missingness difference
plink --bfile after_gender_heterozyg_pihat_mind --exclude missing_snps_1E4.txt --make-bed --out after_gender_heterozyg_pihat_mind_missing1
# create a unique, sorted list of excluded variants
sort -u missing_snps_1E4.txt > VARIANT_TEST_MISSING_SNPS.txt
```
NB: 29 SNPs were filtered out in both cases, but some of the excluded SNPs differ depending on which output files (GCTA or plink) from above were used. 

### 6. Filter Based on Missingness by Haplotype
- This section focuses on variant quality control (QC) based on haplotype missingness. It identifies variants that exhibit significant levels of missingness within certain haplotypes, which may indicate genotyping errors or data quality issues.

- Here, the ```after_gender_heterozyg_pihat_mind_missing1```file genertaed above is used.

```bash
# With GCTA Output
cd ..
mkdir missingness_by_haplotype && cd missingness_by_haplotype
# calculate missingness rates within different haplotype groups
plink --bfile ../missingness_per_variant/after_gender_heterozyg_pihat_mind_missing1 --test-mishap --out missing_hap
# extract SNP IDs with a haplotype-based missingness P-value ≤ 1E-4 (0.0001) 
awk '{if ($8 <= 0.0001) print $9}' missing_hap.missing.hap > missing_haps_1E4.txt
# replace `|` with a newline, listing each SNP ID on its own line (they are listed in pairs)
sed 's/|/\                      
/g' missing_haps_1E4.txt > missing_haps_final.txt 
# create a unique list of the mssingness
sort -u missing_haps_final.txt > HAPLOTYPE_TEST_MISSING_SNPS.txt
# finally, remove SNPs listed in missing_haps_1E4_final.txt due to their significant haplotype-based missingness
plink --bfile ../missingness_per_variant/after_gender_heterozyg_pihat_mind_missing1 --exclude missing_haps_1E4_final.txt --make-bed --out after_gender_heterozyg_pihat_mind_missing2

# declutter the working directoy by saving these to a new dir
mkdir using_plink_output_from_relatedness
mv * missingness_per_variant/using_plink_output_from_relatedness/
```
NB: 15,680 SNPs were removed. 

```bash
# With Plink Output
cd ..
mkdir missingness_by_haplotype && cd missingness_by_haplotype
mkdir using_output_from_plink
# calculate missingness rates within different haplotype groups
plink --bfile ../missingness_per_variant/using_plink_output_from_relatedness/after_gender_heterozyg_pihat_mind_missing1 --test-mishap --out using_output_from_plink/missing_hap
cd using_output_from_plink/
# extract SNP IDs with a haplotype-based missingness P-value ≤ 1E-4 (0.0001) 
awk '{if ($8 <= 0.0001) print $9}' missing_hap.missing.hap > missing_haps_1E4.txt
# replace `|` with a newline, listing each SNP ID on its own line (they are listed in pairs)
sed 's/|/\                      
/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt 
# create a unique list of the mssingness
sort -u missing_haps_1E4_final.txt > HAPLOTYPE_TEST_MISSING_SNPS.txt
# finally, remove SNPs listed in missing_haps_1E4_final.txt due to their significant haplotype-based missingness
plink --bfile ../../missingness_per_variant/after_gender_heterozyg_pihat_mind_missing1 --exclude missing_haps_1E4_final.txt --make-bed --out after_gender_heterozyg_pihat_mind_missing2
```
NB: In this case, 15,602 SNPs were removed instead.

### 8. Hardy-Weinberg Equilibrium (HWE)
In this final QC step, we perform Hardy-Weinberg Equilibrium (HWE) filtering on the remaining SNPs by filtering out variants that deviate significantly from expected allele frequencies, particularly in controls.

```bash
# Since we're using the downstream results obtained after filtering for relatedness with GCTA, only one method is shown here.  
cd ..
mkdir hwe && cd hwe
# restrict the HWE filtering to controls only, which is common in case-control studies to avoid excluding SNPs due to disease-related allele frequency shifts in cases
plink --bfile ../missingness_by_haplotye/after_gender_heterozyg_pihat_mind_missing2 --filter-controls --hwe 1E-4 --write-snplist --out HWE_snps
# extract only the SNPs that passed the HWE filter and create a new binary dataset with these SNPs
plink --bfile ../missingness_by_haplotye/after_gender_heterozyg_pihat_mind_missing2 --extract HWE_snps.snplist --make-bed --out after_gender_heterozyg_pihat_mind_missing3
# Rename the Final Filtered Files for downasstream analysis
mv after_gender_heterozyg_pihat_mind_missing3.bim FILTERED.test.bim
mv after_gender_heterozyg_pihat_mind_missing3.bed FILTERED.test.bed
mv after_gender_heterozyg_pihat_mind_missing3.fam FILTERED.test.fam
mv after_gender_heterozyg_pihat_mind_missing3.log FILTERED.test.log
mv after_gender_heterozyg_pihat_mind_missing3.hh FILTERED.test.hh
```
NB: 916 SNPs from the control dataset were removed here. 

## Imputation
The pipeline below walks through the data-quality checks and formattinmg done pre- and post-imputation of the QC'ed genotype data from above. 

The data is formatted for imputation on the Michigan impuation server

For this:
- ```The HRC-1000G-check-bim-v4.pl``` perl script developed by William Rayner is downloaded from [here](https://www.chg.ox.ac.uk/~wrayner/tools/). The script checks how many variants in the genotype data differ by more than 20% allele frequency.  
- The HRC Panel Haplotype Reference Consortium database for imputation (```HRC.r1 1.GRCh37.wgs.mac5.sites.tab```) is downloaded from [here](http://www.haplotype%20referenceconsortium.org/site).
- A ```filter_variants_maf_r2_updated.R``` R script for fgiltering variants based on MAF, and
- The CheckVCF.py python script is also used (optional) for final QC, provided [here](https://github.com/zhanxw/checkVCF). The reference genome (hs37d5.fa) in FASTA and its index file (hs37d5.fa.fai) must be in the directory when the script is run.

We start by making a new directory and copying the QC'ed data into it.
```bash
# Standing in the main working directory, make new dirs 
mkdir 1.1_post_QC_formatting && cd 1.2_Post_imputation_QC
```

### Post-QC Formatting

```bash
cd 1.1_post_QC_formatting
cp ../01_QC/hwe/FILTERED.test.* .

# Calculate the genotype frequency and run the William Rayner script
FILENAME=FILTERED.test
plink --bfile $FILENAME --freq --out $FILENAME

perl HRC-1000G-check-bim.pl -b $FILENAME.bim -f $FILENAME.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

# NB: This is computationally demanding and can take hours when run on the terminal. 
# Use the SLURM job script below, and it will complete in few minutes.
```

```bash
sbatch perl.sh
```

This outputs the following files, which includes the updates that must be applied to the genotype data:

- ```Run-plink.sh```
- LOG-FILTERED.test-HRC.txt
- Force-Allele1-FILTERED.test-HRC.txt
- Strand-Flip-FILTERED.test-HRC.txt
- ID-FILTERED.test-HRC.txt
- Position-FILTERED.test-HRC.txt
- Chromosome-FILTERED.test-HRC.txt
- Exclude-FILTERED.test-HRC.txt
- FreqPlot-FILTERED.test-HRC.txt

```bash
# Run the generated perl script as follows to apply the corrections:
sh Run-plink.sh
```

This then produces:
- ```FILTERED.test-updated-chr``` plink binary and log files for all 22 chromosomes. 

### Convert to VCF, sort, and compress the outputs

```bash
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
> do
> plink --bfile $FILENAME-updated-chr$chnum --recode vcf --chr $chnum --out $FILENAME$chnum
> done
```

For imputation in the MiCHIGAN imputation Server (MIS), BGEN compressed files are also required.
- We use BCFTOOLS 1.18 for this compression

```bash
module spider bcftools
ml load GCC/12.3.0 BCFtools/1.18

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
> do
> bcftools sort $FILENAME$chnum.vcf -Oz -o pre_impute_$FILENAME$chnum.vcf.gz
> done
```

Finally, do an optional check with the ```checkVCF.py```python script. 
```bash
# This is optional and will depend on the availability of python2.7 version specified by the developer of the script to work
ml spider python
ml GCC/10.2.0 Python/2.7.18 

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
> do
> python2.7 checkVCF.py -r hs37d5.fa -o $FILENAME$chnum pre_impute_$FILENAME$chnum.vcf.gz
done
```
The checkVCF script validates and quality-checks the VCF (Variant Call Format) files to ensure they are suitable for downstream analysis, particularly for the imputation task on MIS. The script performs checkf for:
- Reference Allele Consistency
- Monomorphic Sites (etects sites where there is no variation)
- Indels (identifies non-SNPs)
- Allele Frequency
- Invalid Genotypes


Finally, the files are ready for imputation on [MIS](https://imputationserver.sph.umich.edu/#!). There is a specific guide on the website for this. The results are downloaded after to ```1.2_post_imputation_QC```.

### Post Imputation QC
After imputation, usually not all imputed sites are equal. So, arbitrary cutoff points (softcall and hardcall) for both minor allele frequency and Rsq is used for further QC.

- Softcall: rsq > 0.3 (minimum acceptable level)
- Hardcall: rsq > 0.8 (more stringent)

We also determine what sites pass the MAF and Rsq cutoffs using the R script ```filter_variants_maf_r2_updated.R```provided by the GP2 Parkinson's Disease Genetics course. Note that this script has been modified for this use case.

But, first, the .info files are decompressed
```bash
cd ../1.2_post_imputation_QC
gunzip *.info.gz
# using MAF cutoff of 0.001 and Rsq cutoff of 0.8
Rscript filter_variants_maf_r2_updated.R 0.001 0.8
```

Next, generate the soft and hard call plink binary files using the soft and hard call cutoffs using the bash loop below:

```bash
mkdir soft_call_binaries
mkdir hard_call_binaries
```

```bash
# use "--qual-threshold" to 0.3 to generate softcall

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
> do
> plink --vcf chr$chnum.dose.vcf.gz --qual-scores chr$chnum.info 7 1 1 --qual-threshold 0.3 --make-bed --out soft_call_binaries/chr$chnum  --double-id
> done 

# hard call - change to 0.8
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
> do
> plink --vcf chr$chnum.dose.vcf.gz --qual-scores chr$chnum.info 7 1 1 --qual-threshold 0.8 --make-bed --out hard_call_binaries/chr$chnum  --double-id
> done 
```

```bash
cd soft_call_binaries/
ls | grep ".bim" > merge_list.txt
sed -i'.original' -e 's/.bim//g' merge_list.txt
plink --merge-list merge_list.txt --make-bed --out IMPUTED_SOFTCALLS_MPBC
rm chr*.bim chr*.bed chr*.fam chr*.log chr*.nosex

# Also do this in the hard calls folder
cd ../hard_call_binaries/
ls | grep ".bim" > merge_list.txt
sed -i'.original' -e 's/.bim//g' merge_list.txt
plink --merge-list merge_list.txt --make-bed --out IMPUTED_HARDCALLS_MPBC
rm chr*.bim chr*.bed chr*.fam chr*.log chr*.nosex
```

NB: This project used the softcall
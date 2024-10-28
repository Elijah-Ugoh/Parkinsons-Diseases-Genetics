# Quality Control
This is the first step in the data analysis pipeline, before the PRS is calculated

```bash
mkdir processed_analysis_file/01_QC

# check if plink is installed and load it
module spider plink
module spider PLINK/2.00a3.7
module load GCC/11.3.0  OpenMPI/4.1.4 # load dependencies
module load PLINK/2.00a3.7 # then loda plink
```

## Calculate Genotyping Call Rate per Sample
```bash
plink --bfile MPBC --missing --out processed_analysis_files/01_QC/call_rates
```

### Remove Call Rate Outliers
- Based on the results from the ```processed_analysis_files/1_QC/call_rates.imiss``` file, samples that are not genotyped at the 95% of the tested variance are removed using the ```--mind 0.05``` flag.

```bash
plink --bfile MPBC --mind 0.05 --make-bed --out processed_analysis_files/01_QC/MPBC_call_rates
```

NB: No individual was removed due to missing genotype data 

## Next, Prune the Data for Heterozygosity and Linkage Disequilibrium
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

## Proceed to Sample-level and Gender QC 
```bash
plink --bfile after_het_call_rate --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
```

NB: No gender check problem found. So, the gender obtained from the clinical data correctly aligns with the genotyping data. 
- Make binary files to proceed
```bash
plink --bfile after_het_call_rate --make-bed --out after_gender
```

## Do Ancestry QC
- The samples in the plink binary files is clustered with a reference panel, i.e., [HAPMAP reference populations](http://www.hapmap.org/). This way, all samples in the study that are not from the ancestry of interest (Europe) are filtered out. 
- This comparison with HAPMAP is based on the number of SNPs in the input dataset that overlap with HAPMAP. 

```bash
# move up one directory and make a new directory
cd .. && mkdir ancestry

plink --bfile prunning/after_gender --bmerge ../../HAPMAP_files/HAPMAP_hg19_new --out ancestry/hapmap3_bin_snplis --make-bed
plink --bfile prunning/after_gender --flip ancestry/hapmap3_bin_snplis-merge.missnp --make-bed --out ancestry/after_gender3
plink --bfile ancestry/after_gender3 --bmerge ../../HAPMAP_files/HAPMAP_hg19_new --out ancestry/hapmap3_bin_snplis --make-bed
plink --bfile ancestry/after_gender3 --exclude ancestry/hapmap3_bin_snplis-merge.missnp --out ancestry/after_gender4 --make-bed
plink --bfile ancestry/after_gender4 --bmerge ../../HAPMAP_files/HAPMAP_hg19_new --out ancestry/hapmap3_bin_snplis --make-bed
plink --bfile ancestry/hapmap3_bin_snplis --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out ancestry/pruning
plink --bfile ancestry/hapmap3_bin_snplis --extract ancestry/pruning.prune.in --make-bed --out ancestry/pruned_data
plink --bfile ancestry/pruned_data --het --out ancestry/prunedHet
plink --bfile ancestry/pruned_data --geno 0.01 --out ancestry/pca --make-bed --pca 4
```

### Calculate PCA
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
- "0" for the MPBC population, "1" for the european population, "2" and "3", for the asian and african population, respectively
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

NB: As expected, the MPBC cohort clustered around the population with European ancestry.
 
Finally, on acestry filtering, we keep the samples that cluster with Europeans in the HAPMAP data and make genotype binary files from those.

```bash
plink --bfile ../../prunning/after_gender --keep PCA_filtered_europeans.txt --make-bed --out after_gender_heterozyg_hapmap

# save all the outliers to a different file
cat PCA_filtered_asians.txt PCA_filtered_africans.txt PCA_filtered_mixed.txt > hapmap_outliers.txt 
```

## QC - Filter based on Relatedness
```bash
cd ../../ && mkdir relatedness
plink --bfile ancestry/pca_filtering/after_gender_heterozyg_hapmap --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out relatedness/pruning
plink --bfile ancestry/pca_filtering/after_gender_heterozyg_hapmap --extract relatedness/pruning.prune.in --make-bed --out relatedness/pruned_data
plink --bfile relatedness/pruned_data --het --out relatedness/prunedHet
```
## Run GCTA to Calculate Genetic Relatedness Matrix
Because of the memory requirement for this programme, it is run as a job on the cluster.

Job script:
```bash
#!/bin/bash
#
#SBATCH -J GCTA_GRM
#SBATCH --partition=sens
#SBATCH -o GRM_%j.out
#SBATCH -e GRM_%j.err
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# write to standard output 
cat $0

# load gcta module
module load GCC/12.3.0
module load GCTA/1.94.1

# copy input dta to node local disk
cp -p pruned_data.bim pruned_data.bed pruned_data.fam gcta64 $NAISS_TMP

# change to the execution directory
cd $NAISS_TMP

# run gcta to compute genetic relatedness matrix
# partition the GRM into 3 parts for resource efficiency
gcta64 --bfile pruned_data --make-grm-part 3 1 --out GRM --autosome --maf 0.05 --thread-num 5
gcta64 --bfile pruned_data --make-grm-part 3 2 --out GRM --autosome --maf 0.05 --thread-num 5
gcta64 --bfile pruned_data --make-grm-part 3 3 --out GRM --autosome --maf 0.05 --thread-num 5 

# merge all parts together 
cat GRM.part_3_*.grm.id > GRM.grm.id
cat GRM.part_3_*.grm.bin > GRM.grm.bin
cat GRM.part_3_*.grm.N.bin > GRM.grm.N.bin

# run the relatedness threshold cutoff
gcta64 --grm-cutoff 0.125 --grm GRM --out GRM_0125 --make-grm --thread-num 5

# rescue results to the submission directory
cp -p GRM* $SLURM_SUBMIT_DIR
```





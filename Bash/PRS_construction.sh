#!/bin/bash
### PRS construction Bash script 
### By: Dina Mansour Aly (PhD) Genomics of complex traits
### Karolinska Institute, Sweden


### Input files: snplist.txt (space-delimited file with SNP ids; same in two columns) and PRS_snplist (tab-delimited file with three columns; SNP, Risk-allele, Effect_size).
### Software: PLINK_v1.9

for chr in $(seq 1 22) ; do 

### Step 1: SNP extraction:
../Tools/PLINK_1.9/plink --bfile ../GENOTYPEdata/Binary_genotype_files_for_all_UKBB/Binary_chr"$chr"_ukb --extract snplist.txt --make-bed --out Binary_PRS_chr"$chr"

### Step 2: SNP prunning:
../Tools/PLINK_1.9/plink --bfile ../Binary_PRS_chr"$chr" --indep-pairwise 250 0.5 --out Prunned_PRS_snps
../Tools/PLINK_1.9/plink --bfile ../Binary_PRS_chr"$chr" --extract Prunned_PRS_snps.prune.in --recode tab --out Prunned_PRS_PED_chr"$chr"

### Step 3: PRS score calculation:
../Tools/PLINK_1.9/plink --file ../Prunned_PRS_PED_chr"$chr" --score PRS_snplist.txt

done

#### This outputs the plink.profile that has the raw scores for each individual.
#######################################################################################################################################################

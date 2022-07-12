#! bin/bash
#ssh to server and enter password blindly

cd project3
mkdir Tolani
cd Tolani
mkdir datasets
cd datasets

#download datasets using wget command
wget https://github.com/HackBio-Internship/public_datasets/blob/main/Asia_HLA_Distribution/binary_plink_file/asia.bed.gz?raw=true/ -O asia.bed.gz
wget https://github.com/HackBio-Internship/public_datasets/blob/main/Asia_HLA_Distribution/binary_plink_file/asia.bim?raw=true/ -O asia.bim
wget https://github.com/HackBio-Internship/public_datasets/raw/main/Asia_HLA_Distribution/binary_plink_file/asia.fam
wget https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Global_genome_structure_project/complete_1000_genomes_sample_list_.tsv



#generate eigenvalues with plink
plink --bfile asia --pca 

#perform linkage disequilibrium 

#create a LD pruned set of markers
plink --bfile asia --indep-pairwise 1000 10 0.01 --out prune1

#calculate identity by descent score on the pruned marker list
# DNA segments that are IBD are IBS per definition
plink --bfile asia --extract prune1.prune.in --genome --out ibs1

#cluster individuals into homogeneous groups and perform a multidimensional scaling analysis 
plink --bfile asia --read-genome ibs1.genome --cluster --ppc 1e-3 --cc --mds-plot 2 --out strat1 

#extract SNPs on chromosome 6 and save into a bed file
plink --bfile asia --chr 6  --from-kb 231 --to-kb 171031 --make-bed --out asia_c6

#report pairwise Linkage Disequilibrium (r-squared) for SNPs in this region
plink --bfile asia_c6 --r2 --out asia_c6

#calculate ld matrix and report SNPs with r2 value of 1
plink --bfile asia_c6 --ld-window-r2 1 --out asia_c6f






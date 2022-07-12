#! /bin/bash
#download datasets from github
wget https://github.com/HackBio-Internship/public_datasets/blob/main/vcfs_and_indexes/china.vcf.gz?raw=true/ -O china.vcf.gz
wget https://github.com/HackBio-Internship/public_datasets/blob/main/vcfs_and_indexes/bangladesh.vcf.gz?raw=true/ -O bangladesh.vcf.gz

#unzip datasets
gunzip china.vcf.gz
gunzip bangladesh.vcf.gz

#index files with bcftools
bcftools index china.vcf.gz
bcftools index bangladesh.vcf.gz

#merge files with bcftools
bcftools merge china.vcf bangladesh.vcf -O merged.vcf

#view details of merged files
ls -lh merged.vcf

#! /bin/sh

mkdir stage1
cd stage1
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa

#Question1
#bash code for counting number of DNA sequence
grep -c "^>" DNA.fa 


#Question 2
# bash code for counting total occurrence of A,G,T,C
grep -o 'A\|T\|G\|C' DNA.fa | sort | uniq -c


#Question3
#Setting up a conda environment
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
chmod +x Miniconda3-py38_4.12.0-Linux-x86_64.sh
./Miniconda3-py38_4.12.0-Linux-x86_64.sh
conda activate base
conda --version


#Installing fastqc, multiqc and fastp
conda install -c bioconda fastqc
conda install -c bioconda multiqc
conda install -c bioconda fastp


#Downloading Alsen, Chara and Drysdale R1 & R2 datasets
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true/  -O Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true/ -O Alsen_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true -O Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true -O Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true/ -O Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true/ -O Drysdale_R2.fastq.gz


#creating a folder 'output'
mkdir output


#implementing fastqc on the datasets
fastqc *.fastq.gz -O output/


#implementing multiqc on the qc datasets
multiqc output
 mv multiqc_report.html output 
 mv multiqc_data output

 
#implementing fastp on datasets
fastp -i Alsen_R1.fastq.gz -o output/Alsen_R1.fastq.gz
fastp -i Alsen_R2.fastq.gz -o output/Alsen_R2.fastq.gz
fastp -i Chara_R1.fastq.gz -o output/Chara_R1.fastq.gz
fastp -i Chara_R2.fastq.gz -o output/Chara_R2.fastq.gz
fastp -i Drysdale_R1.fastq.gz -o output/Drysdale_R1.fastq.gz
fastp -i Drysdale_R2.fastq.gz -o output/Drysdale_R2.fastq.gz


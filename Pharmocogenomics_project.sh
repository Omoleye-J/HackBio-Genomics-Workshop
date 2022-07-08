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
plink --bfile asia --pca --chr-set 36 no-xy

#download and install Rstudio on PC
#download datasets including plink outputs from the server into PC using WinSCP

# check the current working directory
getwd()


#sett the directory to the folder containing the downloaded datasets
setwd("Users/User/Downloads/Stage3_datasets")


#create PCA plot
#set eigenvec to pca1
#The file name needs to be enclosed within double quotes, the sep command specifies how columns are defined, and header=F, means that the first line of the file is not column names
pca1 <- read.table("plink.eigenvec",sep=" ",header=F)


#install ggplot on R
#load ggplot
library("ggplot2")

# plot using geom_point with default parameters, this takes the x,y from data_input, columns, V3 and V4.
ggplot(data=pca1, aes(V3,V4)) + geom_point()

#create metatdata using the complete 1000 genome list
metadata <- read.table("complete_1000_genomes_sample_list_.tsv",sep="\t",header=TRUE)

#view header of pca1
head(pca1)
#view header of metadata
head(metadata)

#merge metadata and pca1 
merge_data <- merge(x= pca1,y= metadata, by.x = "V2", by.y = "Sample.name", all = F)

#plot and color by population
ggplot(data=merge_data, aes(V3,V4,color = Population.code)) + geom_point()
#save and close RStudio



#switch to HackBio's linux server
#perform linkage disequilibrium 

#create a LD pruned set of markers
plink --bfile asia --indep-pairwise 50 10 0.2 --out prune1 

#calculate identity by descent score on the pruned marker list
# DNA segments that are IBD are IBS per definition
plink --bfile asia --extract prune1.prune.in --genome --out ibs1

#cluster individuals into homogeneous groups and perform a multidimensional scaling analysis 
plink --bfile asia --read-genome ibs1.genome --cluster --ppc 1e-3 --cc --mds-plot 2 --out strat1 

#adjust filtering values and perform pruning again
plink --bfile asia --indep-pairwise 1000 10 0.01 --out prune2

#extract SNPs on chromosome 6 and save into a bed file
plink --bfile asia --chr 6  --from-kb 231 --to-kb 171031 --make-bed --out asia_c6

#report pairwise Linkage Disequilibrium (r-squared) for SNPs in this region
plink --bfile asia_c6 --r2 --out asia_c6

#calculate ld matrix and report SNPs with r2 value of 1
plink --bfile asia_c6 --ld-window-r2 1 --out asia_c6f

#extract SNPs on chromosome 8 and save into a bed file
plink --bfile asia --chr 8  --from-kb 220 --to-kb 146253 --make-bed --out asia_c8

#report pairwise Linkage Disequilibrium (r-squared) for SNPs in this region
plink --bfile asia_c8 --r2 --out asia_c8

#calculate ld matrix and report SNPs with r2 value of 1
plink --bfile asia_c8 --ld-window-r2 1 --out asia_c8f

#switch back to R studio for further data analysis

#assign asia_c6.ld to tolani
tolani <- read.table('asia_c6.ld', header=T)
#extract SNPS from the specific base pair on chromosome 6
tolani_furtherprune6 <- subset(tolani BP_A > 28477797 &  BP_A < 33448354)
#check dimension of data
dim(tolani_furtherprune6)

#plot the densities of r2 values
boxplot(tolani_furtherprune6$R2)
plot(density(tolani_furtherprune6$R2), main = 'Density Plot of R2')

#further prune out r2 values less than 0.99999
tolani_finalprune6 <- subset(tolani_furtherprune6, R2 >0.99999)

#check dimension of the pruned data
dim(tolani_finalprune6)

#remove repititions and check length of unique values
length(unique(tolani_finalprune6$SNP_B))

#assign asia_c8.ld to tolani
tolani <- read.table('asia_c8.ld', header=T)

#check dimension of data
dim(tolani)

#further prune out r2 values less than 0.99999
tolani_finalprune8 <- subset(tolani, R2 >0.99999)

#check dimension of the pruned data
dim(tolani_finalprune8)

#remove repititions and check length of unique values
length(unique(tolani_finalprune8$SNP_B))


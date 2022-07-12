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

# create colorcoded plot, change axes labels and include title
ggplot(data = merged, aes(V3, V4, color = Population.code)) + geom_point() + xlab("Principal Component 1 (PC1)") + ylab("Principal Component 2 (PC2)") + ggtitle("PCA of selected Asian Populations")



## Multidimensional Scaling Analysis (MDS) Plot
# set strat.mds as newmdsdata
newmdsdata <- read.table("strat.mds", header = TRUE)

# merge metadata and newmdsdata
merged_mds2 <- merge(x = newmdsdata, y = metadata, by.x = "FID", by.y = "Sample.name", all = F)

# non-colorcoded plot of newmdsdata
ggplot(data=newmdsdata, aes(C1,C2)) + geom_point()

# colorcoded plot by population codes
ggplot(data=merged_mds2, aes(C1,C2, color = Population.code)) + geom_point() + ggtitle("Multidimensional Scaling Analysis")

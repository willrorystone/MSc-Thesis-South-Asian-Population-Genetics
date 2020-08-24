####### Script for running anf plotting PCA from a VCF file using SNPRelate #######
######## Written by William Stone, adapted from "http://corearray.sourceforge.net/tutorials/SNPRelate/#:~:text=The%20functions%20in%20SNPRelate%20for,dataset%20from%20specified%20SNP%20eigenvectors" ########

# Clear workspace, format plot output, load SNPRelate package, set working directory
rm(list=ls())
par(mfrow = c(1,1))
par(mar = c(4, 4, 4, 4))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
require(SNPRelate)
setwd("/Users/willstone/Documents/IMPERIAL_COLLEGE/Research_Prodge/1_Results/PCA")
list.files("/Users/willstone/Documents/IMPERIAL_COLLEGE/Research_Prodge/1_Results/PCA")

# attach VCF file
vcf.fn <- "all_pops_22_filtered.recode.vcf"

# reformat VCF to gds file format, summarise the gds file to check
vcf_gds <-  snpgdsVCF2GDS(vcf.fn, "vcf2.gds1", method="biallelic.only")
snpgdsSummary(vcf_gds)

# opening the gds file for analysis
genofile1 <- snpgdsOpen(vcf_gds, allow.duplicate = TRUE)

# conduct LD pruning to LD threshold of 0.5
snpset <- snpgdsLDpruning(genofile1, ld.threshold=0.5)

#unlist the snp set to get every snp ID
snpset.id <- unlist(snpset)

# run PCA analysis using gds file and pruned snp set
pca1 <- snpgdsPCA(genofile1,snp.id=snpset.id, num.thread=2)

# retrieve the explianed variance proportions for each principle component
pca1$varprop

# calculate explained variance as percentage of total variance
pc.percent <- pca1$varprop*100
pc.percent
head(round(pc.percent, 2))

#######  my PC variance proportions as read from standard output
# 1.30 0.37 0.35 0.32 0.30 0.30

#generate a dataframe 'tab' for all eigenvectors from PCA
tab <- data.frame(sample.id = pca1$sample.id,
                  EV1 = pca1$eigenvect[,1], # the first eigenvector
                  EV2 = pca1$eigenvect[,2], # the second eigenvector
                  EV3 = pca1$eigenvect[,3], # third eig
                  EV4 = pca1$eigenvect[,4], #fourth
                  EV5 = pca1$eigenvect[,5], #fifth
                  EV6 = pca1$eigenvect[,6], #sixth
                  stringsAsFactors = FALSE)

# write eignevector data into csv file
write.csv(tab, 'eigenvals.csv', col.names=TRUE)

# read in the eigenvector csv after I manually added the corresponding population to each individual ID
# this allows me to plot PCA results by population
tab <- read.csv("eigenvals.csv", header=TRUE)

#create a matrix of pca plots containing every combo of eigenvectors plotted against each other
#I used this to look for those that explained population structure 
lbls <- paste("PC", 1:6, "\n", format(pc.percent[1:6], digits=2), "%", sep="")
pairs(pca1$eigenvect[,1:6], col=tab$POP, labels=lbls)

#having identified eignevectors associated with structure, I plotted them fully
# creating an object used to define each population as a different shaped data point
shapes <- as.numeric(tab$POP)

#Plot of eigenvector 1 vs eigenvector 2
plot(tab$EV1, tab$EV2, col=as.integer(tab$POP),pch=shapes, xlab="PC1", ylab="PC2")
legend("bottomleft", legend=levels(tab$POP), pch=c(1,2,3,4,5,6), col=1:nlevels(tab$POP), cex=0.75)

#plot of eigenvector 1 vs eigenvector 4
plot(tab$EV6, tab$EV5, col=as.integer(tab$POP), xlab="-PC1", ylab="-PC4")
legend("topright", legend=levels(tab$POP), pch="o", col=1:nlevels(tab$POP))










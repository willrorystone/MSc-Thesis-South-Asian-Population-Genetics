####### Script to analyse and plot data from ADMIXTURE runs ######
#Written by William Stone, adapted from "https://gaworkshop.readthedocs.io/en/latest/contents/07_admixture/admixture.html"

#clear workspace, set the working directory and load dplyr package
rm(list=ls())
setwd("~/Documents/IMPERIAL_COLLEGE/Research_Prodge/1_results/ADMIXTURE")
list.files("~/Documents/IMPERIAL_COLLEGE/Research_Prodge/1_results/ADMIXTURE")
require(dplyr)

#formatting the plot output for a plot with five panels
par(mfrow = c(5,1))
par(mar = c(1.5, 3, 1.5, 3))

#function to label population blocks in of admixture plot
barNaming <- function(vec) {
  retVec <- vec
  for(k in 2:length(vec)) {
    if(vec[k-1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}

#loop to plot the ADMIXTURE Q output files for runs K=2 to K=5
par(mar = c(2.5, 5, 1, 1.5))
for (i in c("all_pops_22_pruned.2.Q", "all_pops_22_pruned.3.Q", "all_pops_22_pruned.4.Q", "all_pops_22_pruned.5.Q", "all_pops_22_pruned.6.Q")){
  table <- read.table(i)
  popmap <- read.table("population_map.txt", col.names=c("Ind","POP"))
  list <- 1:594
  popnumbered <- cbind(popmap, list)
  admixturetab <- cbind(table, popnumbered, by="POP")
  admixturetab <- arrange(admixturetab, desc(V1))
  ordered = admixturetab[order(admixturetab$"POP"),]
  
  barplot(t(as.matrix(ordered)), col=rainbow(6), border=NA,
          names.arg=barNaming(ordered$POP), las=2, cex.names = 0.64)
}

#########
# CV Plot
require("ggplot2")
# read in the csv file containing the CV values I transcribed from standard output after the ADMIXTURE runs.
cv <- read.csv('CV_table.csv', header=TRUE)

#re-format plot output for cv error plot
par(mfrow = c(1,1))
par(mar = c(5, 5, 1, 1.5))

#plot graph of K value agains CV value with a line between the data points
plot <- ggplot(cv, aes(x=K, y=CV))+geom_point()+geom_line()
plot



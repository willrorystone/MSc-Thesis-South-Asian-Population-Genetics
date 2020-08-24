####### Script for analysis of the outputs files from TreeMix ########
#Written by William Stone

#Clear workspace, set working directory and load packages
rm(list=ls())
setwd("~/Documents/IMPERIAL_COLLEGE/Research_Prodge/1_Results/Treemix")
par(mfrow=c(1,1))
require(dplyr)

#accepting code from 'plotting_funcs.R' script supplied by the TreeMix software
source("plotting_funcs.r")

#loop to plot all trees using the output files from treemix runs m=0 to m=6
for (i in c("out_stem", "out_stem1","out_stem2","out_stem3","out_stem4","out_stem5","out_stem6")){
  plot_tree(i)
}

#creating a vector 'f' of the proportions of explained variance associated with each model used
f <- c(get_f("out_stem"),get_f("out_stem1"),get_f("out_stem2"),get_f("out_stem3"),get_f("out_stem4"),get_f("out_stem5"),get_f("out_stem6"))

#putting the proportions of explained variance into a data frame
m <- c(0:6)
data <- data.frame(m,f)
data

#plotting the the proportions of variance in the data explained by each treemix model for model selection
ggplot(data, aes(x=m, y=f))+geom_point()+geom_line()+labs(x="Number of Migration Events")+labs(y="Proportion of Variance Explained")

#plotting the residuals for the m=1 tree that i selected as being the most informative
plot_resid("out_stem1","poporder.txt")



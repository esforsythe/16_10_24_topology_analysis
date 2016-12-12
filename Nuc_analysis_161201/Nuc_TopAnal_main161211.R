#This is the master script for running the topology analysis for many trees. 
#It's designed to be used with the the function, TopAnalFunc, which is compiled in the script, TopAnal_function161107.R

#set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/Nuc_analysis_161201")

#load the needed packages
install.packages("ape")
library(ape)
library(phytools)
library(geiger)
library(lattice)
library(gtable)
install.packages("gridExtra")
library(gridExtra)

###Important note about trees
#For some reason the root function is not working in R.  So instead, trees must be rooted before hand
#I used nw_reroot to root trees before importing them into R
#The grep command for matching the Esal accession is as follows:
#grep -Eoe 'Esal__NC.{20,40}:' RAxML_bipartitions.Atha__NC_000932.1_NP_051* | grep -Eoe 'Esal__NC.*\d' >outgroups

#Update, foound a method that does not require rooting beforehand 

trees<-read.tree("Nuc_cat161211")

Ntrees<-length(trees)

#Plot all the trees
lapply(trees, PlotTreesFunc)

#Now use standard branch lengths
lapply(trees, PlotTreesFunc_branch)

#Apply the function to all the trees
output<-lapply(trees, TopAnalFunc)


#convert the outpt from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=Ntrees, byrow=TRUE))
names(output_df) <- c("Agroup_monophyly", "Crub_Cgrand_monophyly", "Csat_monophyly", "C_group_monophyly", "Topology", "Bootstrap_Support", "Athal_sister")

#For full analysis

#Build a table with the dataframe
pdf("CPtable_161209.pdf", height=30, width=20)
grid.table(output_df)
dev.off()

#make a pie chart of topologies
labels<-names(summary(output_df$Topology))

labels<-paste(labels, summary(output_df$Topology))

pie(summary(output_df$Topology), labels=labels)

#Show the types of non-monophylies
table(output_df$Agroup_monophyly)
table(output_df$Crub_Cgrand_monophyly)
table(output_df$C_group_monophyly)

#look at trees with non-mono A group

#Agroup_nonmono<-subset(output_df, output_df$Agroup_monophyly=="FALSE")
#Athal_sisters_df<-Agroup_nonmono$Athal_sister

###For only the monophyletic trees

#subset the data to only include fully monophyletic A and C groups

subsettedA<-subset(output_df, output_df$Agroup_monophyly=="TRUE")
subsettedCrubCgra<-subset(subsettedA, subsettedA$Crub_Cgrand_monophyly=="TRUE")
subsettedC<-subset(subsettedCrubCgra, subsettedCrubCgra$C_group_monophyly=="TRUE")

#piecharts of monophyletic treees only
labels2<-names(summary(subsettedC$Topology))
labels2<-paste(labels2, summary(subsettedC$Topology))
pie(summary(subsettedC$Topology), labels=labels2)


#Histogram of bootstrap scores
hist(as.numeric(output_df$Bootstrap_Support))

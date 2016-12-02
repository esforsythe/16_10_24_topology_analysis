#This is the master script for running the topology analysis for many trees. 
#It's designed to be used with the the function, TopAnalFunc, which is compiled in the script, TopAnal_function161107.R

#set the working directory
setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/Nuc_analysis_161201")

#load the needed packages
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

trees<-read.tree("Nuctrees_cat")


#Apply the function to all the trees
output<-lapply(trees, TopAnalFunc)


#convert the outpt from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=length(output), byrow=TRUE))
names(output_df) <- c("Agroup_monophyly", "Crub_Cgrand_monophyly", "Csat_monophyly", "C_group_monophyly", "Topology", "Bootstrap_Support")

#For full analysis

#Build a table with the dataframe
pdf("CPtable.pdf", height=100, width=40)
grid.table(output_df)
dev.off()

#make a pie chart of topologies

labels<-names(summary(output_df$Topology))

labels<-paste(labels, summary(output_df$Topology))

pie(summary(output_df$Topology), labels=labels)


#barplot (still need to make into stacked bar chart)
barplot(summary(output_df$Topology))

###For only the monophyletic trees

#subset the data to only include fully monophyletic A and C groups

subsettedA<-subset(output_df, output_df$Agroup_monophyly=="TRUE")
subsettedCrubCgra<-subset(subsettedA, subsettedA$Crub_Cgrand_monophyly=="TRUE")
subsettedCsat<-subset(subsettedCrubCgra, subsettedCrubCgra$Csat_monophyly=="TRUE")
subsettedC<-subset(subsettedCsat, subsettedCsat$C_group_monophyly=="TRUE")

#Show the types of non-monophylies
table(output_df$Agroup_monophyly)
table(output_df$Crub_Cgrand_monophyly)
table(output_df$Csat_monophyly)
table(output_df$C_group_monophyly)

#Make pie chart
labels2<-names(summary(subsettedC$Topology))
labels2<-paste(labels2, summary(subsettedC$Topology))

pie(summary(subsettedC$Topology), labels=labels2)



lapply(trees, PlotTreesFunc)

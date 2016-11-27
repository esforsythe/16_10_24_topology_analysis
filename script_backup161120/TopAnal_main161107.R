#This is the master script for running the topology analysis for many trees. 
#It's designed to be used with the the function, TopAnalFunc, which is compiled in the script, TopAnal_function161107.R

#set the working directory
setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis")

trees<-read.tree("catfile_CP")
list3<-list()

lapply(trees, TopAnalFunc)

list1


#This is the master script for running the topology analysis for many trees. 
#It's designed to be used with the the function, TopAnalFunc, which is compiled in the script, TopAnal_function161107.R

#set the working directory
setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis")

#load the needed packages
library(ape)
library(phytools)
library(geiger)

trees<-read.tree("catfile_CP_edit")

lapply(trees2, PlotTreesFunc)

trees3<-read.tree("catfile_CP_edit")



lapply(trees2, PlotTreesFunc)


trees2<-read.tree("rooted_trees1.newick")

list3<-list()

lapply(trees3, TopAnalFunc2)

#list3



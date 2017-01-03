#This is a script for generating

#set working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/AllNuclear_170102")

#read in dataframe
data<-read.csv("Results_table_170103.csv", header=TRUE)

#filter out "other topologies"
data_noOthers<-subset(data, data$Topology_loose != "Other_topology")

#make pie chart (not including BS scores)
labels_loose<-names(summary(data_noOthers$Topology_loose))
labels_loose<-paste(labels_loose, summary(data_noOthers$Topology_loose))
pie(summary(data_noOthers$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)")

#filter for only well supported topologies (>=70)

data_noOthers70BS<-subset(data_noOthers, Bootstrap_Support >= 70)

#make pie chart (BS scores >= 70)
labels_loose<-names(summary(data_noOthers70BS$Topology_loose))
labels_loose<-paste(labels_loose, summary(data_noOthers70BS$Topology_loose))
pie(summary(data_noOthers70BS$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)\nBS>=70%")




#This is a script for generating figs and calculating D-stat with bootstrap and jackkife
###I call this technique (drumroll......) BLOCK BOOTSTRAPPING

#set working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/AllNuclear_170102")

data<-read.csv("Results_table_170103.csv", header=TRUE)

#filter out "other topologies"
data_noOthers<-subset(data, data$Topology_loose != "Other_topology")

#Find counts of each topology
table<-table(data_noOthers$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#################JACKBOOTING######################
#Split the data into 100 equal blocks
d<-data_noOthers$Topology_loose
split_d<-split(d, ceiling(seq_along(d)/(length(d)/100)))

#Create function for bootstrapping
Boot_Dstat <- function(X=split_d){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(unlist(x.boot))
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 1000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_full<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)


#This is a script for generating

#set working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/AllNuclear_170102")

#Install and load packages
install.packages("boot")
library(boot)

#read in dataframe
data<-read.csv("Results_table_170103.csv", header=TRUE)

#filter out "other topologies"
data_noOthers<-subset(data, data$Topology_loose != "Other_topology")

#make pie chart (not including BS scores)
labels_loose<-names(summary(data_noOthers$Topology_loose))
labels_loose<-paste(labels_loose, summary(data_noOthers$Topology_loose))
pie(summary(data_noOthers$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)")

#Find counts of each topology
table<-table(data_noOthers$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=data_noOthers$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
N<- 10000
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_full<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#############################################################
###Now do the Dstat calculation for individual chromosomes###
#############################################################

#CHROM 1#
chrom1_data<-subset(data_noOthers, Scaffold == 1)

table<-table(chrom1_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom1_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom1<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 2#
chrom2_data<-subset(data_noOthers, Scaffold == 2)

table<-table(chrom2_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom2_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom2<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 3#
chrom3_data<-subset(data_noOthers, Scaffold == 3)

table<-table(chrom3_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom3_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom3<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 4#
chrom4_data<-subset(data_noOthers, Scaffold == 4)

table<-table(chrom4_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom4_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom4<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 5#
chrom5_data<-subset(data_noOthers, Scaffold == 5)

table<-table(chrom1_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom5_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom5<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 6#
chrom6_data<-subset(data_noOthers, Scaffold == 6)

table<-table(chrom6_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom6_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom6<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 7#
chrom7_data<-subset(data_noOthers, Scaffold == 7)

table<-table(chrom7_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom7_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom7<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

#CHROM 8#
chrom8_data<-subset(data_noOthers, Scaffold == 8)

table<-table(chrom8_data$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=chrom8_data$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_chrom8<-c(AC_count, AB_count, Dstat_observed, SD, P_value)

cat_out<-c(out_full, out_chrom1, out_chrom2, out_chrom3, out_chrom4, out_chrom5, out_chrom6, out_chrom7, out_chrom8)


cat_out_df<- data.frame(matrix(unlist(cat_out), nrow=9, byrow=TRUE))
names(cat_out_df) <- c("AC_trees", "AB_trees", "Dstat_obs", "SD", "P_value")

cat_out_df












####filter for only well supported topologies (>=70)

data_noOthers70BS<-subset(data_noOthers, Bootstrap_Support >= 70)

table<-table(data_noOthers70BS$Topology_loose)

AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

Boot_Dstat <- function(X=data_noOthers70BS$Topology_loose){
  x.boot<-sample(X, size=length(X), replace=TRUE)
  table<-table(x.boot)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_temp<-unname(Dstat_temp)
  Dstat_temp
}

N<- 10000
#Perform bootstrap replication
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))

#make pie chart (BS scores >= 70)
labels_loose<-names(summary(data_noOthers70BS$Topology_loose))
labels_loose<-paste(labels_loose, summary(data_noOthers70BS$Topology_loose))
pie(summary(data_noOthers70BS$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)\nBS>=70%")










############failed attempts at bootstrapping############

##Here I tried to replicate the evobiR calcD
#boot.replicate[is.nan(boot.replicate)] <- 0
#StanD<-sd(boot.replicate)
#z<-abs(Dstat_observed/StanD)
#pval<- 2 * (1- pnorm(z))
#pval<-unname(pval)
#pval

###Create D stat function (try1)
#Dstat_func<-function (AB_count, AC_count) {
#  table<-table(data_noOthers$Topology_loose)
#  AB_count<-table["AB_topology"]
#  AC_count<-table["AC_topology"]
#  #stat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
#  #Dstat_temp<-unname(Dstat_temp)
#  return(AB_count, AC_count)
#}

###Create D stat function (try2)
#Dstat_func<-function (tops) {
#  table<-table(tops)
#  AB_count<-table["AB_topology"]
#  AC_count<-table["AC_topology"]
#  Dstat_temp<-(AC_count - AB_count)/(AC_count + AB_count)
#  Dstat_temp<-unname(Dstat_temp)
#  return(Dstat_temp)
#}

###perform bootstrap resample
#boot(data=data_noOthers$Topology_loose, statistic=Dstat_func, R=1000, sim="parametric", stype="i",
#     formula=(AC_count - AB_count)/(AC_count + AB_count))





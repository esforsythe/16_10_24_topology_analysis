#This is a script for generating

#set working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/AllNuclear_170102")

#Install and load packages
install.packages("boot")
library(boot)
install.packages("resample")
library(resample)

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
#Perform bootstrap replication
N<- 10000 #This indicates how many replicates to use
boot.replicate <- replicate(N, Boot_Dstat())
#Calculate P-value
SD<-sd(boot.replicate)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))
out_full<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom1<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom2<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom3<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom4<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom5<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom6<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom7<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

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
out_chrom8<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#Concat the results from all the chromosomes
cat_out<-c(out_full, out_chrom1, out_chrom2, out_chrom3, out_chrom4, out_chrom5, out_chrom6, out_chrom7, out_chrom8)

#convert into dataframe
cat_out_df<- data.frame(matrix(unlist(cat_out), nrow=9, byrow=TRUE))

#append column names and rownames
names(cat_out_df) <- c("AC_trees", "AB_trees", "Dstat_obs", "SD", "Z_score", "P_value")
rownames(cat_out_df) <- c("Full_genome", "Chrom1", "Chrom2", "Chrom3", "Chrom4", "Chrom5", "Chrom6", "Chrom7", "Chrom8")

cat_out_df


#####BLOCK JACKKNIFING######
#Before performing this, you need to make sure 
#the topology data file is sorted by chromosomal 
#position

table<-table(data_noOthers$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

n<-nrow(data_noOthers)
block_size<-round((n/100)-1)
index<-i-1
jackout<-list()

for(i in 1:100) {
 #remove the first block
  if (i == 1){sample_data<-data_noOthers[-c((1):(1+block_size)),]
  }else {
  #remove all the other blocks
    sample_data<-data_noOthers[-c((i*block_size):((i*block_size)+block_size)),]
  }
  #Find the topologies for the rep
  table<-table(sample_data$Topology_loose)
  AB_count<-table["AB_topology"]
  AC_count<-table["AC_topology"]
  BC_count<-table["BC_topology"]
  
  #Calculate the Dstat
  Dstat_rep<-(AC_count - AB_count)/(AC_count + AB_count)
  Dstat_rep<-unname(Dstat_rep)
  #print(Dstat_rep)
  jackout[[i]]<-Dstat_rep
  
}
#clean up the list
jackout<-unlist(jackout)
#Calculate the P-value
SD<-sd(jackout)
z<-(Dstat_observed-0)/SD
P_value<-2*pnorm(-abs(z))

##############################################################################
##############################################################################
##############################################################################
###########filter for only well supported topologies (>=70)###################
##############################################################################
##############################################################################
##############################################################################

data_noOthers70BS<-subset(data_noOthers, Bootstrap_Support >= 70)

#make pie chart (BS scores >= 70)
labels_loose<-names(summary(data_noOthers70BS$Topology_loose))
labels_loose<-paste(labels_loose, summary(data_noOthers70BS$Topology_loose))
pie(summary(data_noOthers70BS$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)\nBS>=70%")


###########D-STAT#############

#Find counts of each topology
table<-table(data_noOthers70BS$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]

#Calculate the Dstat
Dstat_observed<-(AC_count - AB_count)/(AC_count + AB_count)
Dstat_observed<-unname(Dstat_observed)

#Create function for bootstrapping
Boot_Dstat <- function(X=data_noOthers70BS$Topology_loose){
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
out_full<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#############################################################
###Now do the Dstat calculation for individual chromosomes###
#############################################################

#CHROM 1#
chrom1_data<-subset(data_noOthers70BS, Scaffold == 1)

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
out_chrom1<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 2#
chrom2_data<-subset(data_noOthers70BS, Scaffold == 2)

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
out_chrom2<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 3#
chrom3_data<-subset(data_noOthers70BS, Scaffold == 3)

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
out_chrom3<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 4#
chrom4_data<-subset(data_noOthers70BS, Scaffold == 4)

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
out_chrom4<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 5#
chrom5_data<-subset(data_noOthers70BS, Scaffold == 5)

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
out_chrom5<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 6#
chrom6_data<-subset(data_noOthers70BS, Scaffold == 6)

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
out_chrom6<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 7#
chrom7_data<-subset(data_noOthers70BS, Scaffold == 7)

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
out_chrom7<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#CHROM 8#
chrom8_data<-subset(data_noOthers70BS, Scaffold == 8)

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
out_chrom8<-c(AC_count, AB_count, Dstat_observed, SD, z, P_value)

#Concat the results from all the chromosomes
cat_out_df_70BS<-c(out_full, out_chrom1, out_chrom2, out_chrom3, out_chrom4, out_chrom5, out_chrom6, out_chrom7, out_chrom8)

#convert into dataframe
cat_out_df_70BS<- data.frame(matrix(unlist(cat_out), nrow=9, byrow=TRUE))

#append column names and rownames
names(cat_out_df_70BS) <- c("AC_trees", "AB_trees", "Dstat_obs", "SD", "Z_score", "P_value")
rownames(cat_out_df_70BS) <- c("Full_genome", "Chrom1", "Chrom2", "Chrom3", "Chrom4", "Chrom5", "Chrom6", "Chrom7", "Chrom8")

cat_out_df_70BS

################################################################
##########Making tsv files for chrom mapping in coge############
################################################################

data1<-subset(data, Crub_seq != "No_Crub_tips")

data_BC<-subset(data1, Topology_loose == "BC_topology")
data_AC<-subset(data1, Topology_loose == "AC_topology")
data_AB<-subset(data1, Topology_loose == "AB_topology")
data_Othertop<-subset(data1, Topology_loose == "Other_topology")

write.table(data_BC, file='data_BC.tsv', quote=FALSE, sep='\t', row.names=FALSE)
write.table(data_AC, file='data_AC.tsv', quote=FALSE, sep='\t', row.names=FALSE)
write.table(data_AB, file='data_AB.tsv', quote=FALSE, sep='\t', row.names=FALSE)
write.table(data_Othertop, file='data_Othertop.tsv', quote=FALSE, sep='\t', row.names=FALSE)

write.table(data1, file='ALL_genes.tsv', quote=FALSE, sep='\t', row.names=FALSE)


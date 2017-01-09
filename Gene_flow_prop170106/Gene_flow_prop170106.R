#This is a script to calculate the porportion of introgressed loci
#See Green et al. 2010.  There's thorough explanation in the SOM

setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/Gene_flow_prop170106")

install.packages("ape")
library(ape)
library(geiger)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)



trees<-read.tree("All_Nuc_cat170106")
tree<-trees[[6400]]

Drop_tip_top<-function(tree){
if (is.rooted(tree)) {root_tree<-tree} else {
tips<-tree$tip.label
Esal_tip<-grep("Es", tips)
root_tree<-root(tree, Esal_tip, resolve.root=TRUE, edgelabel=TRUE)
}

#store tip names of all relevent species for the droptree
tips<-root_tree$tip.label
Csat_tip<-grep("Cs", tips)
Crub_tip<-grep("Cr", tips)
Cgrand_tip<-grep("Cg", tips)
Athal_tip<-grep("At", tips)
Alyr_tip<-grep("Al", tips)
Bstri_tip<-grep("Bs", tips)
Esal_tip<-grep("Es", tips)

drop_tree<-drop.tip(root_tree, c(Csat_tip, Crub_tip, Cgrand_tip), trim.internal=TRUE, subtree=FALSE)

tips2<-drop_tree$tip.label
Athal_tip<-grep("At", tips2)
Alyr_tip<-grep("Al", tips2)
Bstri_tip<-grep("Bs", tips2)
Esal_tip<-grep("Es", tips2)

#test if all A group seqs are monophyletic and print results
Agroup_mono<-is.monophyletic(phy=drop_tree, c(Athal_tip, Alyr_tip))

#test if all Bstri and Alyr seqs are monophyletic and print results
B_Alyr_mono<-is.monophyletic(phy=drop_tree, c(Bstri_tip, Alyr_tip))

#test if all A group seqs are monophyletic and print results
B_Athal_mono<-is.monophyletic(phy=drop_tree, c(Athal_tip, Bstri_tip))

###########################################################################################
############                          Topology analysis                        ############
###########################################################################################

#Store the correct topology
if(Agroup_mono) {final_topology = "A_group_topology"} else if(B_Alyr_mono) {final_topology = "B_Alyr_topology"} else if (B_Athal_mono) {final_topology = "B_Athal_topology"} else {final_topology = "Other_topology"}

#Store the node representing the MRCA of each potential clade
#This is the node label at which the crucial BS score resides
AA_MRCA<-getMRCA(phy=drop_tree, c(Athal_tip, Alyr_tip))
BAlyr_MRCA<-getMRCA(phy=drop_tree, c(Bstri_tip, Alyr_tip))
BAthal_MRCA<-getMRCA(phy=drop_tree, c(Bstri_tip, Athal_tip))

#retrieve the supporting BS score
if(final_topology == "A_group_topology") {BS_score = (drop_tree$node.label[(AA_MRCA - length(drop_tree$tip.label))])} else if(final_topology == "B_Alyr_topology") {BS_score = (drop_tree$node.label[(BAlyr_MRCA - length(drop_tree$tip.label))])} else if(final_topology == "B_Athal_topology") {BS_score = (drop_tree$node.label[(BAthal_MRCA - length(drop_tree$tip.label))])} else {BS_score = "BS_scoreNA"} 

#Find Athal tip names 
Athal_tip2<-grep("At", tips2, value=TRUE)
if (length(Athal_tip2) == 0)
{Athal_tip_name = "No_Athal_tips"} else if (length(Athal_tip2)==1)
{Athal_tip_name = Athal_tip2} else if (length(Athal_tip2)==2)
{Athal_tip_name = paste(Athal_tip2[1], Athal_tip2[2])} else if (length(Athal_tip2)==3)
{Athal_tip_name = paste(Athal_tip2[1], Athal_tip2[2], Athal_tip2[3])}

#return results
return(c(Athal_tip_name, final_topology, BS_score))
}

#Apply the function to all the trees
output<-lapply(trees, Drop_tip_top)

#Count number of trees 
Ntrees<-length(trees)

#convert the outpt from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=Ntrees, byrow=TRUE))
names(output_df) <- c("Athal_tip", "Topology", "BS_score")

#For full analysis

#export a CSV file
write.csv(output_df, file = "test.csv")


#Build a table with the dataframe
pdf("retdups_table_161230.pdf", height=30, width=20)
grid.table(output_df)
dev.off()

#make a pie chart of topologies
labels<-names(summary(output_df$Topology))
labels<-paste(labels, summary(output_df$Topology))
pie(summary(output_df$Topology), labels=labels, main="Nuclear Topologies")

#filter out "other topologies"
data_noOthers<-subset(output_df, output_df$Topology != "Other_topology")

#Find counts of each topology
table<-table(data_noOthers$Topology)
AA_count<-table["A_group_topology"]
BAlyr_count<-table["B_Alyr_topology"]
BAthal_count<-table["B_Athal_topology"]

#Calculate the denominator for calculating the propotion of introgressed loci
Denom<-AA_count-BAlyr_count
Denom<-unname(Denom)
#Denom = 6152   (17-02-16)

##########  Now I'm going to switch over to a different WD to calculate the numerator  ######
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/AllNuclear_170102")
data<-read.csv("Results_table_170103.csv", header=TRUE)
data_noOthers<-subset(data, data$Topology_loose != "Other_topology")
table<-table(data_noOthers$Topology_loose)
AB_count<-table["AB_topology"]
AC_count<-table["AC_topology"]
BC_count<-table["BC_topology"]
Num4propGF<-(AC_count - AB_count)
Num4propGF<-unname(Num4propGF)

########## Now I'm switching back to the WD inwhich this script lives  ######
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/Gene_flow_prop170106")

#Finally, calculating the proportion of introgressed loci! 
#This should be multiplied by 100 for a percetage
PropIG<-(Num4propGF/Denom)

#Using Athal as the "P3"
#PropIG = 0.03722367    (17-01-06)

#Using Alyr as the "P3"
#PropIG = 0.0371874    (17-01-06)


#This is a script for compiling a topology analysis function

setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis")


#load the needed packages
library(ape)
library(phytools)
library(geiger)

tree<-read.tree("RAxML_bipartitions.Athal__AT5G27950.1.phy")

tree2<-read.tree("RAxML_bipartitions.Athal__AT5G61130.1.phy")

testtree<-read.tree("testtree")


TopAnalFunc<-function(tree){
  

  #list the tips that are Csativa accessions
  #Csat_tips<-grep("Csativa", tips)
  
  #store the names of Csativa tips in excess of one (these will be 'dropped' later)
#  Csat_tips4drop<-if(Csat_counts == 3) {tail(Csat_tips, 2)} else if (Csat_counts == 2) {tail(Csat_tips, 1)}
  
  
  #test if Csativa seqs are monophyletic and print results
  #Csat_mono<-is.monophyletic(phy=tree, c("Csat_tips1"))
 # if(Csat_mono) {print("Csativa monophyletic")} else {print("Csativa non-monophyletic")}
  
  
  #if necessary, drop excess Csativa tips and create new "droptree" for downstream analysis
 # if (Csat_counts > 1) {tree<-drop.tip(tree,Csat_tips4drop)} else {tree<-tree}
  
  
  #store tip names of all relevent species for the droptree
  tips<-testtree$tip.label
  Csat_tips<-grep("Csat", tips)
  Crub_tip<-grep("Crub", tips)
  Cgrand_tip<-grep("Cgra", tips)
  Athal_tip<-grep("Atha", tips)
  Alyr_tip<-grep("Alyr", tips)
  Bstri_tip<-grep("Bstr", tips)
  Esal_tip<-grep("Esal", tips)
  
  
  ###########################################################################################
  ############                    A group monophyly test                         ############
  ###########################################################################################
  
  
  #test if A group seqs are monophyletic and print results
  Agroup_mono<-is.monophyletic(phy=testtree, c(Athal_tip, Alyr_tip))
  if(Agroup_mono) {Agroupmono = "A group monophyletic"} else {Agroupmono = "A group non-monophyletic"}
  
  
  ###########################################################################################
  ############                    Crub-Cgra group monophyly test                         ############
  ###########################################################################################

  
  #test if Crub-Cgrand seqs are monophyletic and print results
  CrubCgrand_mono<-is.monophyletic(phy=testtree, c(Crub_tip, Cgrand_tip))
  if(CrubCgrand_mono) {Capgroupmono = "Crub-Cgrand monophyletic"} else {Capgroupmono = "Crub-Cgrand non-monophyletic"}
  
  ###########################################################################################
  ############                     C.sativa paralog monophyly test                         ############
  ###########################################################################################
  
  #test if Csativa seqs are monophyletic and print results
  Csat_mono<-is.monophyletic(phy=testtree, c("Csat_tips"))
  if(Csat_mono) {Csatmono= "Csativa monophyletic"} else {Csatmono= "Csativa non-monophyletic"}


  ###########################################################################################
  ############                    Full C group monophyly test                         ############
  ###########################################################################################
  
  
  #test if all C group seqs are monophyletic and print results
  Cgroup_mono<-is.monophyletic(phy=testtree, c(Crub_tip, Cgrand_tip, Csat_tips))
  if(Cgroup_mono) {Cgroupmono = "C group monophyletic"} else {Cgroupmono = "Csat_keeper_tip non-monophyletic"}
  
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
  rooted_tree<-root(testtree, Esal_tip, resolve.root = TRUE)

  
  ###########################################################################################
  ############                   Topology analysis of keeper trees               ############
  ###########################################################################################
  
  #Check which clade is monophyletic
  BC_clade<-is.monophyletic(phy=rooted_tree, c(Crub_tip, Cgrand_tip, Csat_tips, Bstri_tip))
  AC_clade<-is.monophyletic(phy=rooted_tree, c(Crub_tip, Cgrand_tip, Csat_tips, Athal_tip, Alyr_tip))
  AB_clade<-is.monophyletic(phy=rooted_tree, c(Athal_tip, Alyr_tip, Bstri_tip))
  
  #Store the correct topology
  if(BC_clade) {final_topology = "BC_topology"} else if(AC_clade) {final_topology = "AC_topology"} else if (AB_clade) {final_topology = "AB_topology"}

  #Store the node representing the MRCA of each potential clade
  #This is the node label at which the crucial BS score resides
  BC_MRCA<-getMRCA(phy=rooted_tree, c(Crub_tip, Cgrand_tip, Csat_tips, Bstri_tip))
  AC_MRCA<-getMRCA(phy=rooted_tree, c(Crub_tip, Cgrand_tip, Csat_tips, Athal_tip, Alyr_tip))
  AB_MRCA<-getMRCA(phy=rooted_tree, c(Athal_tip, Alyr_tip, Bstri_tip))
 
 #plot.phylo(rooted_tree, show.node.label=TRUE)
 
 #retrieve the supporting BS score
 if(BC_clade) {BS_score = (rooted_tree$node.label[(BC_MRCA - length(rooted_tree$tip.label))])} else if(AC_clade) {BS_score = (rooted_tree$node.label[(AC_MRCA - length(rooted_tree$tip.label))])} else if(AB_clade) {BS_score = (rooted_tree$node.label[(AB_MRCA - length(tree$tip.label))])} 

#list3<-list()
#list2<-list()

list3<-c(Agroupmono, Capgroupmono, Csatmono, Cgroupmono, final_topology, BS_score)
#list2$BS_score<-BS_score

}



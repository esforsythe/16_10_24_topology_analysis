setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis")

#load the needed packages
library(ape)
library(phytools)
library(geiger)

#read in specific tree 
tree<-read.tree("RAxML_bipartitions.Athal__AT5G27950.1.phy")


TreePooper<-function(tree){
  
  tips<-tree$tip.label
  #count how many Csativa accessions are present
  Csat_counts<-length(grep("Csativa", tips))
  
  #list the tips that are Csativa accessions
  Csat_tips<-grep("Csativa", tips)
  
  #store the names of Csativa tips in excess of one (these will be 'dropped' later)
  Csat_tips4drop<-if(Csat_counts == 3) {tail(Csat_tips, 2)} else if (Csat_counts == 2) {tail(Csat_tips, 1)}
  
  
  #test if Csativa seqs are monophyletic and print results
  Csat_mono<-is.monophyletic(phy=tree, c("Csat_tips1"))
  if(Csat_mono) {print("Csativa monophyletic")} else {print("Csativa non-monophyletic")}
  
  
  #if necessary, drop excess Csativa tips and create new "droptree" for downstream analysis
  if (Csat_counts > 1) {Csat_dropped_tree<-drop.tip(tree,Csat_tips4drop)} else {Csat_dropped_tree<-tree}
  
  
  #store tip names of all relevent species for the droptree
  tips_droptree<-Csat_dropped_tree$tip.label
  Csat_keeper_tip<-grep("Csativa", tips_droptree)
  Crub_tip<-grep("Crub", tips_droptree)
  Cgrand_tip<-grep("Cgrand", tips_droptree)
  Athal_tip<-grep("Athal", tips_droptree)
  Alyr_tip<-grep("Alyr", tips_droptree)
  Bstri_tip<-grep("Bstri", tips_droptree)
  Esal_tip<-grep("Esal", tips_droptree)
  
  
  ###########################################################################################
  ############                    A group monophyly test                         ############
  ###########################################################################################
  
  
  #test if A group seqs are monophyletic and print results
  Agroup_mono<-is.monophyletic(phy=Csat_dropped_tree, c(Athal_tip, Alyr_tip))
  if(Agroup_mono) {print("A group monophyletic")} else {print("A group non-monophyletic")}
  
  
  ###########################################################################################
  ############                    Crub-Cgra group monophyly test                         ############
  ###########################################################################################
  
  #store tip names
  Crub_tip<-grep("Crub", tips_droptree)
  Cgrand_tip<-grep("Cgrand", tips_droptree)
  
  
  #test if Crub-Cgrand seqs are monophyletic and print results
  CrubCgrand_mono<-is.monophyletic(phy=Csat_dropped_tree, c(Crub_tip, Cgrand_tip))
  if(CrubCgrand_mono) {print("Crub-Cgrand monophyletic")} else {print("Crub-Cgrand non-monophyletic")}
  
  ###########################################################################################
  ############                    Full C group monophyly test                         ############
  ###########################################################################################
  
  
  #test if all C group seqs are monophyletic and print results
  Cgroup_mono<-is.monophyletic(phy=Csat_dropped_tree, c(Crub_tip, Cgrand_tip, Csat_keeper_tip))
  if(CrubCgrand_mono) {print("C group monophyletic")} else {print("Csat_keeper_tip non-monophyletic")}
  
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
  rooted_tree<-root(Csat_dropped_tree, Esal_tip, resolve.root = TRUE)
  is.rooted(rooted_tree)
  
  ###########################################################################################
  ############                   Topology analysis of rooted tree               ############
  ###########################################################################################
  
  #find the node representing the most recent common ancestor of the C clade
  MRCA_C_clade<-getMRCA(tree, c(Csat_keeper_tip, Crub_tip, Cgrand_tip))
  
  #find the node that is sister to the C-clade
  sister_node<-getSisters(tree, MRCA_C_clade)
  
  #find the tips that are in that sister clade
  sister_tips<-tips(tree, sister_node)
  
  #use the paste function to combide all the sister tips into one string so that I can grep that long string
  sistertips_string<-paste(c(sister_tips), collapse=', ')
  
  
  #assess whether sister clade containes Athal or Bstri
  topology<-if(grep("Athal", sistertips_string)) {("A-C topology")} else if (grep("Bstri", sistertips_string)) {("B-C topology")} 
  
  #override the above topology if the sister clad contains both Athal and Bstri.  If the sister clade contains both, it means that the tree topology is the A-B topology
  topology_final<-if (grepl("(?=.*Athal)(?=.*Bstri)", sistertips_string, perl = TRUE)) {("A-B topology")} else {(topology)}
  
  
  #COMMAND FOR GET BOOTSTRAP:
  tree$node.label[(connecting_node - length(tree$tip.label))]
  
  plot(topology_final)
}










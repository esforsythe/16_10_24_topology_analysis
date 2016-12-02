#This is a script for compiling a topology analysis function

setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/CP_analysis_161201")

#This is the simplest function I can think of.  I need to test if my functions are working on many trees
PlotTreesFunc<-function(tree){
  #tips<-tree$tip.label
  #Esal_tip<-grep("Esal", tips)
  #tree<-root(tree, Esal_tip)
  #tree<-compute.brlen(tree, 100)
  plot.phylo(tree)
  nodelabels(tree$node.label)
}

PlotTreesFunc_branch<-function(tree){
  #tips<-tree$tip.label
  #Esal_tip<-grep("Esal", tips)
  #tree<-root(tree, Esal_tip)
  tree<-compute.brlen(tree, 100)
  plot.phylo(tree)
  nodelabels(tree$node.label)
}


#Here I'm gonna try to root all the trees
RootTreesFunc<-function(tree){
  tips<-tree$tip.label
  Esal_tip<-grep("Esal", tips)
  #sink("rooted_trees.newick", append=TRUE)
  rooted<-(root(tree, Esal_tip))
  write.tree(rooted, file = "rooted_trees1.newick", append = TRUE, digits = 5, tree.names = FALSE)
  #sink()
}




#######   BELOW IS FOR TESTING FUNCTION   #######
#making a test tree
#text_tree<-"((Atha__NC_000932.1_NP_051039.1_2:0.00584499745201570806,(Csat__psbA:0.00095886724450046444,Bstri_NP_051039:0.00386671023621103574)38:0.00000044403951708354)33:0.00096083638464386761,(Esal__NC_028170.1_YP_009175655.1_2:0.01092018124633497775,(Cgra__NC_028517.1_YP_009182871.1_2:0.00000044403951708354,Crub__NC_027693.1_YP_009161901.1_2:0.00000044403951708354)97:0.00095549074757341238)47:0.00096953690884424531,Alyr__LN877383.1_CUA65384.1_2:0.00192104722868210532);"
#text_tree<-"(Cgra__NC_028517.1_YP_009182872.1_3:0.00000075476092639459,(((Csat__matK:0.01691398929133055076,(Atha__NC_000932.1_NP_051040.2_3:0.01078436325752899459,Alyr__LN877383.1_CUA65385.1_3:0.00830728909446542228)76:0.00187919888266765978)24:0.00051981043633657223,Bstri_NP_051040:0.01136718501490044213)30:0.00100034132674294680,Esal__NC_028170.1_YP_009175656.1_3:0.03489016167715788819)100:0.01747755328112032824,Crub__NC_027693.1_YP_009161902.1_3:0.00000075476092639459);"


#reading test tree
#tree<-read.tree(text=text_tree)

#plotting test tree
#plot.phylo(tree)

#root test tree
#tips<-tree$tip.label
#Esal_tip<-grep("Esal", tips)
#rooted_tree<-root(tree, Esal_tip, resolve.root = TRUE)

#look at bootstrap support values
#rooted_tree$node.label

#set branch lengths to 1
#rooted_tree<- compute.brlen(rooted_tree, 100)

#plotting test tree
#plot.phylo(rooted_tree)
#nodelabels(rooted_tree$node.label)

#check if tree is bifrucating
#is.binary.tree(tree)



#Here is a command for clearing all object. Be careful with it.
#rm(list=ls())

tree<-read.tree("testtree")

####### END TESTING #######



TopAnalFunc<-function(tree){
  
  #store tip names of all relevent species for the droptree
  tips<-tree$tip.label
  Csat_tip<-grep("Csat", tips)
  Crub_tip<-grep("Crub", tips)
  Cgrand_tip<-grep("Cgra", tips)
  Athal_tip<-grep("Atha", tips)
  Alyr_tip<-grep("Alyr", tips)
  Bstri_tip<-grep("Bstr", tips)
  Esal_tip<-grep("Esal", tips)
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
 #tree<-root(tree, Esal_tip, resolve.root=TRUE)
  
  
  ###########################################################################################
  ############                    A group monophyly test                         ############
  ###########################################################################################
  
  
  #test if A group seqs are monophyletic and print results
  Agroup_mono<-is.monophyletic(phy=tree, c(Athal_tip, Alyr_tip))
  if(Agroup_mono) {Agroupmono = "A group monophyletic"} else {Agroupmono = "A group non-monophyletic"}
  
  
  ###########################################################################################
  ############                    Crub-Cgra group monophyly test                         ############
  ###########################################################################################

  
  #test if Crub-Cgrand seqs are monophyletic and print results
  CrubCgrand_mono<-is.monophyletic(phy=tree, c(Crub_tip, Cgrand_tip))
  if(CrubCgrand_mono) {Capgroupmono = "Crub-Cgrand monophyletic"} else {Capgroupmono = "Crub-Cgrand non-monophyletic"}
  
  ###########################################################################################
  ############                     C.sativa paralog monophyly test                         ############
  ###########################################################################################
  
  #test if Csativa seqs are monophyletic and print results
  Csat_mono<-is.monophyletic(phy=tree, c("Csat_tips"))
  if(Csat_mono) {Csatmono= "Csativa monophyletic"} else {Csatmono= "Csativa non-monophyletic"}


  ###########################################################################################
  ############                    Full C group monophyly test                         ############
  ###########################################################################################
  
  
  #test if all C group seqs are monophyletic and print results
  Cgroup_mono<-is.monophyletic(phy=tree, c(Crub_tip, Cgrand_tip, Csat_tip))
  if(Cgroup_mono) {Cgroupmono = "C group monophyletic"} else {Cgroupmono = "C group non-monophyletic"}
  
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
  #rooted_tree<-root(tree, Esal_tip, resolve.root = TRUE)

  
  ###########################################################################################
  ############                   Topology analysis of keeper trees               ############
  ###########################################################################################
  
  #Check which clade is monophyletic
  BC_clade<-is.monophyletic(phy=tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))
  AC_clade<-is.monophyletic(phy=tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))
  AB_clade<-is.monophyletic(phy=tree, c(Athal_tip, Alyr_tip, Bstri_tip))
  
  #Store the correct topology
  if(BC_clade) {final_topology = "BC_topology"} else if(AC_clade) {final_topology = "AC_topology"} else if (AB_clade) {final_topology = "AB_topology"} else {final_topology = "Other_topology"}

  #Store the node representing the MRCA of each potential clade
  #This is the node label at which the crucial BS score resides
  BC_MRCA<-getMRCA(phy=tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))
  AC_MRCA<-getMRCA(phy=tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))
  AB_MRCA<-getMRCA(phy=tree, c(Athal_tip, Alyr_tip, Bstri_tip))
 
 #plot.phylo(rooted_tree, show.node.label=TRUE)
 
 #retrieve the supporting BS score
 if(BC_clade) {BS_score = (tree$node.label[(BC_MRCA - length(tree$tip.label))])} else if(AC_clade) {BS_score = (tree$node.label[(AC_MRCA - length(tree$tip.label))])} else if(AB_clade) {BS_score = (tree$node.label[(AB_MRCA - length(tree$tip.label))])} else {BS_score = "BS_scoreNA"} 

 #Investigating the non-monophyletic topologies
 Athal_sister<-c(tips(tree, getSisters(tree, Athal_tip, mode="number")))
 
 
 return(c(Agroup_mono, CrubCgrand_mono, Cgroup_mono, final_topology, BS_score))
          #Athal_sister))
}



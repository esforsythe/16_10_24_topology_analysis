#This is a script for compiling a topology analysis function

setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/Nuc_analysis_161201")

#This is the simplest function I can think of.  I need to test if my functions are working on many trees
PlotTreesFunc<-function(tree){
  tree<-compute.brlen(tree, 100)
  tips<-tree$tip.label
  Esal_tip<-grep("Esal", tips)
  root_tree<-root(tree, Esal_tip, resolve.root=TRUE)
  root_tree$node.label=tree$node.label
  plot.phylo(root_tree, cex=0.7)
  nodelabels(root_tree$node.label, cex=0.7)
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

tree<-read.tree(
  text="((Csativa__XP_010457102:0.01518881265333598893,(Csativa__XP_010480051:0.00662865210348756104,Csativa__XP_010474770:0.00671291887992642774)52:0.00324457624170141744)100:0.02191034610995688109,((Crub__Carubv10010070m:0.00000153396898798239,Cgrand__Cagra.1968s0103:0.00395412686650943156)100:0.03470315157896354530,(Bstri__Bostr.5325s0084:0.03216375026092795769,(Alyr__16054647:0.02326242938724031040,Athal__AT1G01630:0.02159315289028859372)100:0.04204974080296691036)39:0.00542206353861695092)50:0.00982388583721963138,Esal__Thhalv10008535m:0.13360687503277937638);
")

####### END TESTING #######

TopAnalFunc<-function(tree){
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
  #normalize the branch lengths because short branches are problematic
  #tree<-compute.brlen(tree, 100)
  
  tips<-tree$tip.label
  Esal_tip<-grep("Esal", tips)
  root_tree<-root(tree, Esal_tip, resolve.root=TRUE)
  
  #After rooting, the BS scores are all messed up. Below, I take the scores from the unrooted tree and put them on the rooted tree.  I need to double-check to make sure this is working correctly! 
  
  root_tree$node.label=tree$node.label
  #root_tree$node.label
  
  
  #store tip names of all relevent species for the droptree
  tips2<-root_tree$tip.label
  Csat_tip<-grep("Csat", tips2)
  Crub_tip<-grep("Crub", tips2)
  Cgrand_tip<-grep("Cgra", tips2)
  Athal_tip<-grep("Atha", tips2)
  Alyr_tip<-grep("Alyr", tips2)
  Bstri_tip<-grep("Bstr", tips2)
  Esal_tip<-grep("Esal", tips2)
  
  #Collapse poorly supported branches
  #This is a function that someone wrote.  I found it at:
  #http://stackoverflow.com/questions/34403957/how-to-collapse-branches-in-a-phylogenetic-tree-by-the-label-in-their-nodes-or-l
 #not working right now
 # root_tree<-di2multi4node(root_tree)
  
  ###########################################################################################
  ############                    A group monophyly test                         ############
  ###########################################################################################
  
  
  #test if A group seqs are monophyletic and print results
  Agroup_mono<-is.monophyletic(phy=root_tree, c(Athal_tip, Alyr_tip))
  if(Agroup_mono) {Agroupmono = "A group monophyletic"} else {Agroupmono = "A group non-monophyletic"}
  
  
  ###########################################################################################
  ############                    Crub-Cgra group monophyly test                         ############
  ###########################################################################################
  
  
  #test if Crub-Cgrand seqs are monophyletic and print results
  CrubCgrand_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip))
  if(CrubCgrand_mono) {Capgroupmono = "Crub-Cgrand monophyletic"} else {Capgroupmono = "Crub-Cgrand non-monophyletic"}
  
  ###########################################################################################
  ############                     C.sativa paralog monophyly test                         ############
  ###########################################################################################
  
  #test if Csativa seqs are monophyletic and print results
  Csat_mono<-is.monophyletic(phy=root_tree, c(Csat_tip))
  if(Csat_mono) {Csatmono= "Csativa_monophyletic"} else {Csatmono= "Csativa_non_monophyletic"}
  
  
  ###########################################################################################
  ############                    Full C group monophyly test                         ############
  ###########################################################################################
  
  
  #test if all C group seqs are monophyletic and print results
  Cgroup_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip))
  if(Cgroup_mono) {Cgroupmono = "C group monophyletic"} else {Cgroupmono = "C group non-monophyletic"}
  
  
  ###########################################################################################
  ############                   Topology analysis of keeper trees               ############
  ###########################################################################################
  
  #Check which clade is monophyletic
 
  BC_clade<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))

  AC_clade<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))

  AB_clade<-is.monophyletic(phy=root_tree, c(Athal_tip, Alyr_tip, Bstri_tip))
  
  #Store the correct topology (STRICT: requires that A, B, and/or C clade are perfectly monogomous)
  if(BC_clade & Cgroup_mono) {final_topology = "BC_topology"} else if(AC_clade & Cgroup_mono & Agroup_mono) {final_topology = "AC_topology"} else if (AB_clade & Agroup_mono) {final_topology = "AB_topology"} else {final_topology = "Other_topology"}
  
  #Store the correct topology (LOOSE: Topology analysis with less-strict monophyly requirements)
  if(BC_clade) {final_topology_loose = "BC_topology"} else if(AC_clade) {final_topology_loose = "AC_topology"} else if (AB_clade) {final_topology_loose = "AB_topology"} else {final_topology_loose = "Other_topology"}
  
  #Store the node representing the MRCA of each potential clade
  #This is the node label at which the crucial BS score resides
  BC_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Bstri_tip))
  AC_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip, Athal_tip, Alyr_tip))
  AB_MRCA<-getMRCA(phy=root_tree, c(Athal_tip, Alyr_tip, Bstri_tip))
  
  #plot.phylo(root_tree, show.node.label=TRUE)
  
  #retrieve the supporting BS score
  if(final_topology_loose == "BC_topology") {BS_score = (root_tree$node.label[(BC_MRCA - length(root_tree$tip.label))])} else if(final_topology_loose == "AC_topology") {BS_score = (root_tree$node.label[(AC_MRCA - length(root_tree$tip.label))])} else if(final_topology_loose == "AB_topology") {BS_score = (root_tree$node.label[(AB_MRCA - length(root_tree$tip.label))])} else {BS_score = "BS_scoreNA"} 
  
  ###Investigating the non-monophyletic topologies
  #What is sister to A. thaliana?
  Athal_sisters<-c(tips(root_tree, getSisters(root_tree, Athal_tip, mode="number")))
  Athal_sis_count<-length(Athal_sisters)
  if(Athal_sis_count==1) 
  {Athal_sis_paste = Athal_sisters} else if(Athal_sis_count==2) 
  {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2])} else if(Athal_sis_count==3) 
  {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3])} else if(Athal_sis_count==4) 
  {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4])} else if(Athal_sis_count==5) 
  {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4], Athal_sisters[5])} else {Athal_sis_paste = "Many_sisters"}
 
  #What is sister to the C.grand - C. rubella clade?
  CrubCgra_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip))
  CrubCgra_sisters<-c(tips(root_tree, getSisters(root_tree, CrubCgra_MRCA, mode="number")))
  CrubCgra_sis_count<-length(CrubCgra_sisters)
  if(CrubCgra_sis_count==1) 
  {CrubCgra_sis_paste = CrubCgra_sisters} else if(CrubCgra_sis_count==2) 
  {CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2])} else if(CrubCgra_sis_count==3) 
  {CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3])} else if(CrubCgra_sis_count==4) 
  {CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3], CrubCgra_sisters[4])} else if(CrubCgra_sis_count==5) 
  {CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3], CrubCgra_sisters[4], CrubCgra_sisters[5])} else {CrubCgra_sis_paste = "Many_sisters"}
  
 #plot.phylo(root_tree, show.node.label=TRUE)
 
 
  
  return(c(Agroup_mono, CrubCgrand_mono, Csat_mono, Cgroup_mono, final_topology, final_topology_loose, BS_score, Athal_sis_paste, CrubCgra_sis_paste))
}



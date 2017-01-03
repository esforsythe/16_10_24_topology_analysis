#This is a script for compiling a topology analysis function

setwd("~/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/CP_analysis_161201")

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
text="(Esal__NC_028170.1_YP_009175659.1_6:0.00000042022745740336,(Cgra__NC_028517.1_YP_009182875.1_6:0.00000042022745740336,(Atha__NC_000932.1_NP_051043.1_6:0.00000042022745740336,((Bstri_NP_051043:0.02875295565858896479,Csat__psbI:0.00933938042193883412)30:0.00000042022745740336,Alyr__LN877383.1_CUA65388.1_6:0.00000042022745740336)10:0.00000042022745740336)0:0.00000042022745740336)10:0.00000042022745740336,Crub__NC_027693.1_YP_009161905.1_6:0.00000042022745740336);
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
  root_tree<-root(tree, Esal_tip, resolve.root=TRUE, edgelabel=TRUE)
  
  #After rooting, the BS scores are all messed up. Below, I take the scores from the unrooted tree and put them on the rooted tree.  I need to double-check to make sure this is working correctly! 
  
 # root_tree$node.label=tree$node.label
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
  Csat_mono<-is.monophyletic(phy=root_tree, c("Csat_tips"))
  if(Csat_mono) {Csatmono= "Csativa monophyletic"} else {Csatmono= "Csativa non-monophyletic"}


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
  
 #Investigating the non-monophyletic topologies
 Athal_sisters<-c(tips(root_tree, getSisters(root_tree, Athal_tip, mode="number")))
 Athal_sis_count<-length(Athal_sisters)
 if(Athal_sis_count==1) 
    {Athal_sis_paste = Athal_sisters} else if(Athal_sis_count==2) 
        {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2])} else if(Athal_sis_count==3) 
            {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3])} else if(Athal_sis_count==4) 
                {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4])} else if(Athal_sis_count==5) 
                    {Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4], Athal_sisters[5])} else {Athal_sis_paste = Many_sisters}
 
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
 
 
 
 #find the name of the Athal tip(s)
 Athal_tip2<-grep("At", tips2, value=TRUE)
 if (length(Athal_tip2) == 0)
 {Athal_tip_name = "No_Athal_tips"} else if (length(Athal_tip2)==1)
 {Athal_tip_name = Athal_tip2} else if (length(Athal_tip2)==2)
 {Athal_tip_name = paste(Athal_tip2[1], Athal_tip2[2])} else if (length(Athal_tip2)==3)
 {Athal_tip_name = paste(Athal_tip2[1], Athal_tip2[2], Athal_tip2[3])}
 
 #find the name of the Crub tip(s)
 Crub_tip2<-grep("Cr", tips2, value=TRUE)
 if (length(Crub_tip2) == 0)
 {Crub_tip_name = "No_Crub_tips"} else if (length(Crub_tip2)==1)
 {Crub_tip_name = Crub_tip2} else if (length(Crub_tip2)==2)
 {Crub_tip_name = paste(Crub_tip2[1], Crub_tip2[2])} else if (length(Crub_tip2)==3)
 {Crub_tip_name = paste(Crub_tip2[1], Crub_tip2[2], Crub_tip2[3])}
 
 
 return(c(Athal_tip_name, Crub_tip_name, Agroup_mono, CrubCgrand_mono, Csat_mono, Cgroup_mono, final_topology, final_topology_loose, BS_score))
 #Athal_sis_paste, CrubCgra_sis_paste))
}



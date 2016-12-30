
#set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/RetDups_161213")


#load the tree
#tree<-read.tree("RAxML_bipartitions.197.fasta.phy")

#bad tree
tree<-read.tree(
  text="((Al16056296:0.02020789127187572892,(AtAT1G50370:0.01020848708191410964,((Cs010500597:0.01072540316942311583,(Cs010479490:0.00531405032717378357,Cs010461883:0.00505037496084351570)59:0.00175709866842810943)93:0.00740131147735283244,((Cr10009861:0.00000122416941185337,Cr10009924:0.06004453337246020567)41:0.00002760583436008093,Cg0860s0009:0.00222765298919147103)100:0.01477114112670410533)78:0.00607033763205656364)30:0.00267391425505425972)54:0.00496257291335305394,(Es10012080:0.02014233310303580396,(Es10021218:0.03970824946514577708,((Al16035856:0.00858423422270451417,AtAT3G19980:0.00768625385329832315)92:0.01116422061867688650,((Cs010487957:0.00589729473554996206,Cs010466166:0.00210628990960991117)100:0.01773038773075499841,(Cg0316s0069:0.02018260974294596236,Bs19424s0665:0.00611192337277987729)88:0.00480026963228620367)96:0.00927698303012416264)81:0.01246785650453173348)100:0.06134530904120149425)67:0.01084841293945711625,Bs26675s0138:0.01030456751422858912);")

#good tree
tree<-read.tree(
  text="(Cr10008569:0.00468282090420213837,((Bs15697s0399:0.01677484614195268869,(((Es10016381:0.05398522695737057592,((AtAT2G34940:0.03974505475831915241,Al16056355:0.02887199802043730404)94:0.01341205980929299942,(((Cg0993s0005:0.00437417468304125961,Cr10024993:0.00493702189921782844)100:0.04038141989710786889,Cs010509684:0.02200666803589884646)96:0.00867038179455673942,Bs23794s0609:0.02040302005503837057)100:0.02174489728399313371)100:0.01855269516638761301)100:0.19464646726415377187,Es10007059:0.02898551807730217686)99:0.02338780154194186264,(Al16038880:0.01651057157828074440,AtAT1G30900:0.03008566000603494009)91:0.00946897471965313375)54:0.00336938384550336794)90:0.00694056664682361508,(Cs010478532:0.01063949917459242056,Cs010460907:0.01377500543503866297)100:0.01799081422679884440)100:0.01885820403267945977,Cg1508s0107:0.01006451231069125513);
")

#BS debug tree 
tree<-read.tree(
  text="(Al16046805:0.02151999879561814627,((Es10007402:0.05149637094984380048,(Es10028585:0.01854350125185300416,((AtAT4G02290:0.02652812623090433641,Al16047732:0.01452782676346934913)92:0.00497778618478454411,((((Cs010456052:0.00064572451763527294,Cs010456051:0.00000195809054684368)100:0.00527513805751742698,(Cs010422611:0.00000195809054684368,Cs010422612:0.00000195809054684368)100:0.00387544078492336758)100:0.01510840443387256697,(Cr10000741:0.00533068260121769347,Cg0753s0009:0.00384008143586313721)100:0.02484910618776361990)100:0.02079819868951569492,Bs16335s0022:0.01571388746784372140)99:0.01101894025011419809)100:0.03511076597586191422)100:0.17064354211939899986)99:0.02832484043025221276,(Bs4824s0048:0.01920704727219900779,(Cs010481246:0.03228961428521013449,(Cr10011914:0.00271303524365368873,Cg0566s0020:0.00131461156353825559)100:0.02494710055904493518)85:0.00460279412807546029)95:0.00862302908823404414)99:0.01201443765595408379,AtAT1G02800:0.01770422986556156181);
")

#BS debug tree TWO 
tree<-read.tree(
  text="(AtAT1G03440:0.04187114736661239067,((((Cs010482019:0.00892043219658301273,Cs010457352:0.01435909806981883544)100:0.03064189833869035429,(Cg0897s0036:0.00948660195299155795,Cr10009374:0.00843410616892512066)100:0.04479526052592704144)100:0.01799473060997335352,Bs26675s0409:0.02800045298680505640)95:0.01359664144612090717,(Es10007833:0.05296089315262035579,(Es10028712:0.07209245987680516776,((Al16052825:0.02856084135001413227,AtAT4G03010:0.03853733051962800887)98:0.01418603820720442006,(Bs2983s0141:0.02700137720334141317,(Cs010430633:0.04686687234906011801,(Cg2552s0026:0.00960537314393852641,Cr10003741:0.00260572075228526117)100:0.03463138745472209984)89:0.01408569171677174700)65:0.00736980102632710677)99:0.04028433347754449556)100:0.14153114841339650698)77:0.01424026568817608751)95:0.01769571435218073030,Al16055687:0.02174917917597851691);
")

#new bad tree
tree<-read.tree(
  text="(AtAT4G14900:0.04008419396965915565,(((Cs010449973:0.03665365337269542445,(Cr10004542:0.00443075038548971633,Cg2641s0005:0.00758902659326777043)100:0.02998424201653140941)89:0.01034982607738572230,Bs0597s0119:0.02450593630240286563)45:0.00633034450822827100,(Es10024868:0.08270945915136457272,(Es10020440:0.22163269662277917949,((Cg1189s0022:0.00630882942844495850,Cr10015242:0.00330432597010198951)100:0.02385457879406777115,(Cs010511354:0.04016760734812741679,(Bs19424s0409:0.03291457933612540465,(Al16051682:0.01941687644435566296,AtAT3G22440:0.01934887297824286478)97:0.01369438471327411380)33:0.00555438774823359144)61:0.01140108755778668909)99:0.03529516997949427165)100:0.13199523602852064608)92:0.03067898637436113421)100:0.01520343529155917736,Al16040497:0.02992910645218251758);")

#new new bad tree
tree<-read.tree(
  text="(((Cg424:0.00083794383615743448,Cr408:0.00184797329647181643)100:0.02406522073183675267,((Cs1006:0.00000132114869771561,Cs1003:0.00000132114869771561)100:0.01176170245164121636,Cs1004:0.01211493163422812380)99:0.00524486435021751570)87:0.00452230956219556220,((((Es414:0.05275161243444976417,((Bs408:0.01878629976476501109,(Cs1002:0.03666316697503224642,(Cg423:0.00345904104551136371,Cr407:0.00187147136699476312)100:0.02569209198044747344)100:0.01807002220204762324)84:0.00608522077439483788,At461:0.03403276578671626806)99:0.01173501891893004821)100:0.05445483946438144790,(Al420:0.03715381575081574267,Cs1005:0.05158681800660384897)100:0.33558204061195134882)100:0.06755086646527501404,Es413:0.04475580169867881403)96:0.02157639733576732430,(Al419:0.00680102797085907311,(At462:0.00000132114869771561,At463:0.00000132114869771561)100:0.01869691827183433269)100:0.01051170892634918641)60:0.00357239585372689551,Bs407:0.01891340623498140652);")


PlotTreesFunc<-function(tree){
  plot.phylo(tree, cex=0.7, show.node.label=TRUE)
}


###Fuction for splitting trees into sub trees
Split_trees<-function(tree){
  
#midpoint root the tree
root_tree<-midpoint(tree, node.labels="support")

#root_tree$node.label=tree$node.label
#normalize branch lengths
#this is required to avoid the "incorrect number of dimensions" error from treeSlice
root_tree<-compute.brlen(root_tree, 100)

#Slice tree at the root (actually very very close to the root)
#This outputs a multiphylo of the two sub trees
sliced_trees<-treeSlice(root_tree, 0.01, trivial=FALSE, prompt=FALSE)

#This if/else loop is to verify that the split worked and yeilded two sub trees
#if it didn't yeild two subtrees, that means that one of the accessions came out sister to all others.  These trees should be discarded. 
if (length(sliced_trees) == 2){
alltips<-root_tree$tip.label
tips1<-sliced_trees[[1]]$tip.label
tips2<-sliced_trees[[2]]$tip.label

subtree1<-drop.tip(root_tree, setdiff(alltips, tips2))
subtree2<-drop.tip(root_tree, setdiff(alltips, tips1))

#print trees
write.tree(c(subtree1, subtree2), append=TRUE)
}
#^^^ end of if loop ^^^

#vvv end of split function vvv
}

#example tree
tree<-read.tree(
  text="(Es10012080:100,(Bs26675s0138:100,((((Cg0860s0009:100,(Cr10009924:100,Cr10009861:100)30:100)78:100,((Cs010461883:100,Cs010479490:100)NA:100,Cs010500597:100)54:100)93:100,AtAT1G50370:100)59:100,Al16056296:100)100:100)41:100)92:99.99;")

#missing BS score tree
tree<-read.tree(
  text="(Es10001934:100,(((Cs010417438:100,Cs010472674:100)100:100,(Cg2858s0006:100,Cr10024513:100)71:100)70:100,(Bs8111s0008:100,(AtAT2G26450:100,Al16056254:100)NA:100):100)100:100)100:99.99;")


TopAnalFunc<-function(tree){
  
  ###########################################################################################
  ############                    Root the tree by the Esal tip                  ############
  ###########################################################################################
  
  #normalize the branch lengths because short branches are problematic
  #tree<-compute.brlen(tree, 100)
  
  #tips<-tree$tip.label
  #Esal_tip<-grep("Es", tips)
  
  #If there is no Esal tip, midpoint root the subtree
  #if(Esal_tip == 0) {root_tree<-midpoint(tree)} else { 
  
  #root_tree<-root(tree, Esal_tip, resolve.root=TRUE)}
  
  
  #After rooting, the BS scores are all messed up. Below, I take the scores from the unrooted tree and put them on the rooted tree.  I need to double-check to make sure this is working correctly! 
  
  #root_tree$node.label=tree$node.label
  #root_tree$node.label
  
  root_tree<-tree

  #store tip names of all relevent species for the droptree
  tips2<-root_tree$tip.label
  Csat_tip<-grep("Cs", tips2)
  Crub_tip<-grep("Cr", tips2)
  Cgrand_tip<-grep("Cg", tips2)
  Athal_tip<-grep("At", tips2)
  Alyr_tip<-grep("Al", tips2)
  Bstri_tip<-grep("Bs", tips2)
  Esal_tip<-grep("Es", tips2)
  
  #Collapse poorly supported branches
  #This is a function that someone wrote.  I found it at:
  #http://stackoverflow.com/questions/34403957/how-to-collapse-branches-in-a-phylogenetic-tree-by-the-label-in-their-nodes-or-l
  #not working right now
  # root_tree<-di2multi4node(root_tree)
  
  ###########################################################################################
  ############                    A group monophyly test                         ############
  ###########################################################################################
  
  
  #test if A group seqs are monophyletic and print results
  #if(Athal_tip != 0 & Alyr_tip != 0){
  Agroup_mono<-is.monophyletic(phy=root_tree, c(Athal_tip, Alyr_tip))
  #if(Agroup_mono) {Agroupmono = "A group monophyletic"} else {Agroupmono = "A group non-monophyletic"}
#} else {Agroup_mono == "NA"}
  
  ###########################################################################################
  ############                    Crub-Cgra group monophyly test                         ############
  ###########################################################################################
  
  
  #test if Crub-Cgrand seqs are monophyletic and print results
#if(Crub_tip != 0 & Cgrand_tip != 0){
  CrubCgrand_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip))
 # if(CrubCgrand_mono) {Capgroupmono = "Crub-Cgrand monophyletic"} else {Capgroupmono = "Crub-Cgrand non-monophyletic"}
#} else {CrubCgrand_mono == "NA"}

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
#if(Crub_tip != 0 & Cgrand_tip != 0 & Csat_tip != 0){
  Cgroup_mono<-is.monophyletic(phy=root_tree, c(Crub_tip, Cgrand_tip, Csat_tip))
 # if(Cgroup_mono) {Cgroupmono = "C group monophyletic"} else {Cgroupmono = "C group non-monophyletic"}
#} else {Cgroup_mono == "NA"}
  
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
  
  #.phylo(root_tree, show.node.label=TRUE)
  
  #retrieve the supporting BS score
  if(final_topology_loose == "BC_topology") {BS_score = (root_tree$node.label[(BC_MRCA - length(root_tree$tip.label))])} else if(final_topology_loose == "AC_topology") {BS_score = (root_tree$node.label[(AC_MRCA - length(root_tree$tip.label))])} else if(final_topology_loose == "AB_topology") {BS_score = (root_tree$node.label[(AB_MRCA - length(root_tree$tip.label))])} else {BS_score = "BS_scoreNA"} 
  
  ###Investigating the non-monophyletic topologies
  ##What is sister to A. thaliana?
  #Athal_sisters<-c(tips(root_tree, getSisters(root_tree, Athal_tip, mode="number")))
  #Athal_sis_count<-length(Athal_sisters)
  #if(Athal_sis_count==1) 
  #{Athal_sis_paste = Athal_sisters} else if(Athal_sis_count==2) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2])} else if(Athal_sis_count==3) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3])} else if(Athal_sis_count==4) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4])} else if(Athal_sis_count==5) 
  #{Athal_sis_paste = paste(Athal_sisters[1], Athal_sisters[2], Athal_sisters[3], Athal_sisters[4], Athal_sisters[5])} else {Athal_sis_paste = "Many_sisters"}
  
  ###What is sister to the C.grand - C. rubella clade?
  #CrubCgra_MRCA<-getMRCA(phy=root_tree, c(Crub_tip, Cgrand_tip))
  #CrubCgra_sisters<-c(tips(root_tree, getSisters(root_tree, CrubCgra_MRCA, mode="number")))
  #CrubCgra_sis_count<-length(CrubCgra_sisters)
  #if(CrubCgra_sis_count==1) 
  #{CrubCgra_sis_paste = CrubCgra_sisters} else if(CrubCgra_sis_count==2) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2])} else if(CrubCgra_sis_count==3) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3])} else if(CrubCgra_sis_count==4) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3], CrubCgra_sisters[4])} else if(CrubCgra_sis_count==5) 
  #{CrubCgra_sis_paste = paste(CrubCgra_sisters[1], CrubCgra_sisters[2], CrubCgra_sisters[3], CrubCgra_sisters[4], CrubCgra_sisters[5])} else {CrubCgra_sis_paste = "Many_sisters"}
  
  #.phylo(root_tree, show.node.label=TRUE)
  
  
  
  return(c(Agroup_mono, CrubCgrand_mono, Csat_mono, Cgroup_mono, final_topology, final_topology_loose, BS_score
           #, Athal_sis_paste, CrubCgra_sis_paste
           ))
}




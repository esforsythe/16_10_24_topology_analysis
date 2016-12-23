#set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/RetDups_161213")


trees<-read.tree("cat_retdups_161223")

out_trees_string<-sapply(trees, Split_trees)

out_trees<-read.tree(text=c(out_trees_string))


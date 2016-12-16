
#set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/RetDups_161213")

#load some packages
install.packages("ape")
library(ape)
install.packages("phytools")
library(phytools)
library(geiger)
library(phangorn)

#load the tree
tree<-read.tree("RAxML_bipartitions.197.fasta.phy")


#Function for midpoint rooting
#obtained from http://grokbase.com/t/r/r-sig-phylo/109268tgx8/midpoint-rooting
midpoint <- function(tree){
  dm = cophenetic(tree)
  tree = unroot(tree)
  rn = max(tree$edge)+1
  maxdm = max(dm)
  ind = which(dm==maxdm,arr=TRUE)[1,]
  tmproot = Ancestors(tree, ind[1], "parent")
  tree = phangorn:::reroot(tree, tmproot)
  edge = tree$edge
  el = tree$edge.length
  children = tree$edge[,2]
  left = match(ind[1], children)
  tmp = Ancestors(tree, ind[2], "all")
  tmp= c(ind[2], tmp[-length(tmp)])
  right = match(tmp, children)
  if(el[left]>= (maxdm/2)){
    edge = rbind(edge, c(rn, ind[1]))
    edge[left,2] = rn
    el[left] = el[left] - (maxdm/2)
    el = c(el, maxdm/2)
  }
  else{
    sel = cumsum(el[right])
    i = which(sel>(maxdm/2))[1]
    edge = rbind(edge, c(rn, tmp[i]))
    edge[right[i],2] = rn
    eltmp = sel[i] - (maxdm/2)
    # el = c(el, sel[i] - (maxdm/2))
    el = c(el, el[right[i]] - eltmp)
    el[right[i]] = eltmp
  }
  tree$edge.length = el
  tree$edge=edge
  tree$Nnode = tree$Nnode+1
  phangorn:::reorderPruning(phangorn:::reroot(tree, rn))
}

#midpoint root the tree
root_tree<-midpoint(tree)

#Slice tree at the root (or very very close to the root)
#This outputs a multiphylo of the two sub trees
sliced_trees<-treeSlice(root_tree, 0.01, trivial=TRUE, prompt=FALSE)

#TO DO:
#figure out how to make this into a fuction an lapply it to many trees
#Need to figure out have to concatinate all the multiphylos into a large multiphylo


#Also, the BS are messed up.  This is a related problem to the one I faced with the single copy gene families.

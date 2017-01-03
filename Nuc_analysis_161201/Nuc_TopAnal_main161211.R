#This is the master script for running the topology analysis for many trees. 
#It's designed to be used with the the function, TopAnalFunc, which is compiled in the script, TopAnal_function161107.R

#set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/Brassicaceae_Phylo/16_10_24_topology_analysis/Nuc_analysis_161201")

#load the needed packages
install.packages("ape")
library(ape)
library(geiger)
library(lattice)
library(gtable)
install.packages("gridExtra")
library(gridExtra)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)

###Important note about trees
#For some reason the root function is not working in R.  So instead, trees must be rooted before hand
#I used nw_reroot to root trees before importing them into R
#The grep command for matching the Esal accession is as follows:
#grep -Eoe 'Esal__NC.{20,40}:' RAxML_bipartitions.Atha__NC_000932.1_NP_051* | grep -Eoe 'Esal__NC.*\d' >outgroups

#Update, foound a method that does not require rooting beforehand 

trees<-read.tree("Nuc_cat161211")

Ntrees<-length(trees)

#Plot all the trees
lapply(trees, PlotTreesFunc)

#Now use standard branch lengths
lapply(trees, PlotTreesFunc_branch)

#Apply the function to all the trees
output<-lapply(trees, TopAnalFunc)


#convert the outpt from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=Ntrees, byrow=TRUE))
names(output_df) <- c("Athl_seq", "Crub_seq", "Agroup_monophyly", "Crub_Cgrand_monophyly", "Csat_monophyly", "C_group_monophyly", "Topology", "Topology_loose", "Bootstrap_Support")

#For full analysis

#export a CSV file
write.csv(output_df, file = "Nuc_topAnal170102.csv")

#Build a table with the dataframe
pdf("Nuc_table_161209.pdf", height=30, width=20)
grid.table(output_df)
dev.off()

#make a pie chart of topologies (STRICT)
labels<-names(summary(output_df$Topology))

labels<-paste(labels, summary(output_df$Topology))

pie(summary(output_df$Topology), labels=labels, main="Nuclear Topologies (strict)")

#make a pie chart of topologies (LOOSE)
labels_loose<-names(summary(output_df$Topology_loose))

labels_loose<-paste(labels_loose, summary(output_df$Topology_loose))

pie(summary(output_df$Topology_loose), labels=labels_loose, main="Nuclear Topologies (loose)")


#Show the types of non-monophylies
table(output_df$Agroup_monophyly)
table(output_df$Crub_Cgrand_monophyly)
table(output_df$Csat_monophyly)
table(output_df$C_group_monophyly)

######Venn Diagram###########
install.packages('VennDiagram')
library(VennDiagram)

#find the rows that have each group non-monophyletic (aka FALSE)
A_false<-nrow(subset(output_df, output_df$Agroup_monophyly == FALSE))
CrubCgra_false<-nrow(subset(output_df, output_df$Crub_Cgrand_monophyly == FALSE))
Csat_false<-nrow(subset(output_df, output_df$Csat_monophyly == FALSE))
C_false<-nrow(subset(output_df, output_df$C_group_monophyly == FALSE))

#Find how the falses overlap eachother
A_C_false<-nrow(subset(output_df, output_df$Agroup_monophyly == FALSE & output_df$C_group_monophyly == FALSE ))
A_Csat_false<-nrow(subset(output_df, output_df$Agroup_monophyly == FALSE & output_df$Csat_monophyly == FALSE))
C_Csat_false<-nrow(subset(output_df, output_df$C_group_monophyly == FALSE & output_df$Csat_monophyly == FALSE))
A_C_Csat_false<-nrow(subset(output_df, output_df$Agroup_monophyly == FALSE & output_df$C_group_monophyly == FALSE & output_df$Csat_monophyly == FALSE))

#Plot the venn diagram
#Note: grid.newpage was throwing an error so I ran dev.off() and it totally fixed it! Go figure! 
grid.newpage()
draw.triple.venn(area1 = A_false, area2 = C_false, area3 = Csat_false,
                 n12 = A_C_false, n23 = C_Csat_false, n13 = A_Csat_false, n123 = A_C_Csat_false,
                    category=c("Agroup", "Cgroup", "Csat_paralogs"),
                        lty = "blank", 
                            fill = c("skyblue", "pink1", "mediumorchid"))
                                                                                                                                                                 fill = c("skyblue", "pink1", "mediumorchid"))))


#look at trees with non-mono A group

#Agroup_nonmono<-subset(output_df, output_df$Agroup_monophyly=="FALSE")
#Athal_sisters_df<-Agroup_nonmono$Athal_sister

###For only the monophyletic trees

#subset the data to only include fully monophyletic A and C groups

subsettedA<-subset(output_df, output_df$Agroup_monophyly=="TRUE")
subsettedCrubCgra<-subset(subsettedA, subsettedA$Crub_Cgrand_monophyly=="TRUE")
subsettedC<-subset(subsettedCrubCgra, subsettedCrubCgra$C_group_monophyly=="TRUE")

#piecharts of monophyletic treees only
labels2<-names(summary(subsettedC$Topology))
labels2<-paste(labels2, summary(subsettedC$Topology))
pie(summary(subsettedC$Topology), labels=labels2)


###Histogram of bootstrap scores
#filter for the ones that aren't NA
boots_only<-subset(output_df, output_df$Bootstrap_Support != "BS_scoreNA")

#plot  the histogram
full_hist<-hist(as.numeric(as.vector(boots_only$Bootstrap_Support)), breaks=10)

#Split the topologies
boots_AC<-subset(boots_only, boots_only$Topology_loose == "AC_topology")
AC_hist<-hist(as.numeric(as.vector(boots_AC$Bootstrap_Support)), breaks=10)

boots_BC<-subset(boots_only, boots_only$Topology_loose == "BC_topology")
BC_hist<-hist(as.numeric(as.vector(boots_BC$Bootstrap_Support)), breaks=10)

boots_AB<-subset(boots_only, boots_only$Topology_loose == "AB_topology")
AB_hist<-hist(as.numeric(as.vector(boots_AB$Bootstrap_Support)), breaks=10)

plot(full_hist, main="distibution of BS scores", xlab= "BS score", border="black", xlim=c(0,100))
plot(AC_hist, main="distibution of BS scores (AC topology)", xlab= "BS score", border="red", xlim=c(0,100))
plot(BC_hist, main="distibution of BS scores (BC topology)", xlab= "BS score", border="blue", xlim=c(0,100))
plot(AB_hist, main="distibution of BS scores (AB topology)", xlab= "BS score", border="green", xlim=c(0,100))



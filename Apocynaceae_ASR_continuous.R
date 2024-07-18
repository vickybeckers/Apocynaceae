# Packages
library(phytools)
library(geiger)
library(viridis) # colour palette for the colour blind
library(tidyverse) # Used for changing colours in plot
library(dplyr) # needed for function n_distinct
library(ape) # needed for read.tree


# Data
setwd("C:/Users/Vicky Beckers")

tree<-read.tree=("Name.tre") # molecular phylogeny from Fishbein et al., 2018
#treeL<-ladderize(tree, right=FALSE) # ladderize if preferred

data  <- read.table("20220729_WA_dataset.txt", row.names=1, header=T)


# Data formatting
data <- as.data.frame(data$Trait, stringAsFactors=FALSE, row.names=trait$Species_name) # 147 observations
dataS <- subset(data, data[,1]!='NA') # # remove rows with missing data, 1 = column number in data
td <- treedata(treeL, dataS)

treeTD<-td$phy # here is the tree
dataTD<-td$data # here is the data


# Trait mapping
traitC<-setNames(dataTD[,1],rownames(dataTD))
obj<-contMap(treeTD,traitC,plot=FALSE)


# Visualisation
n_cols<-n_distinct(traitC)
contmap_obj_viridis<-setMap(obj, viridis(n_cols, direction=-1)) # colour scale
plot(contmap_obj_viridis, legend=0.7*max(nodeHeights(treeTD)), fsize=c(0.3,0.6), lwd=c(2,3), outline=FALSE)
# Adjusted code from Github page of Joyce Chery which is used for Pace et al., 2022 (https://doi.org/10.1007/s40415-021-00764-2)


# Packages
library(phytools)
library(geiger)
library(viridis) # colour palettes for the colour blind
library(ape) # needed for read.tree


# Data
setwd("C:/Users/Vicky Beckers")

tree<-read.tree=("Name.tre") # molecular phylogeny from Fishbein et al., 2018
is.ultrametric(tree)
#treeU<-force.ultrametric(tree, method="extend") # for when is ultrametric but r doesn't read it as such
#treeULT<-ladderize(tree, right=FALSE)

data  <- read.table("20220729_WA_dataset.txt", row.names=1, header=T) 


# Data formatting
LA <- as.data.frame(cbind(rownames(data), data$Trait), row.names = FALSE)
LAS <- subset(LA, LA[,2]!='NA') # remove rows with missing data, 2 = column number in data

species <- LAS$V1 # define species
LAR <- LAS$V2 # define trait
names(LAR) <- species

td <- treedata(treeULT, LAR); # prune tree to fit trait dataset
treeTD<-td$phy # define where tree is


# Ancestral state reconstruction of discrete traits
er  <- fitMk(treeTD, LAR, model="ER", ordered=TRUE, pi="estimated")
ard <- fitMk(treeTD, LAR, model="ARD", ordered=TRUE, pi="estimated")

  ## chi-square test
1 - pchisq(2*abs(ard$logLik- er$logLik), 1) # if p is < .05 select ARD model

  ## stochastic mapping
simmap.trees <- make.simmap(treeTD, LAR, model = "ARD", pi = "estimated", nsim=1000)

  ## set colours depending on number of states
cols<-c("#fde725", "#440154")
names(cols)<-c("Absent", "Present")

  ## visualisation (after blog Revell http://blog.phytools.org/2022/06/how-to-plot-tip-node-labels-without.html)
dataRW<-as.factor(data$Trait)
activity.pattern<-setNames(dataRW,rownames(data))
head(activity.pattern)
plotTree(treeTD,ftype="i",fsize=0.3,lwd=1,offset=1, type="fan")
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
points(pp$xx,pp$yy)
h<-max(nodeHeights(treeTD))
plotTree(treeTD,ftype="i",fsize=0.2,lwd=1,offset=5, type="fan")
points(pp$xx[1:Ntip(treeTD)]+0.01*h,pp$yy[1:Ntip(treeTD)], pch=19,
       col=cols[activity.pattern[treeTD$tip.label]], cex=1.7)
for(i in 1:treeTD$Nnode){
  plotrix::floating.pie(pp$xx[i+Ntip(treeTD)],
                        pp$yy[i+Ntip(treeTD)],radius=0.025*h,
                        x=obj$ace[i,],col=cols, pch=16, border="transparent")
}
legend("bottomleft",levels(activity.pattern),pch=16,col=cols,
       pt.cex=2,bty="n")
# Packages
library(phytools)
library(ape) # needed for read.tree


# Data
setwd("C:/Users/Vicky Beckers")
tree<-read.tree=("Name.tre") # molecular phylogeny from Fishbein et al., 2018
data <- read.table("20230207_WA_Morpho_dataset.txt", row.names=1, header=T)


# Data formatting
df1 <- data[c(9,11,38,39)] # subsetting useful columns
df2 <- na.omit(df1) # delete rows with missing data
tree2 <- treedata(tree, df2) # prune tree to fit data table 
tree3 <- tree2$phy ## here is the tree pruned


# Loop with Pagel's test for correlated evolution for entire dataset
df.res <- data.frame()
for(i in 1:ncol(df2)){
  # reformat character 1
  datum <- as.data.frame(cbind(rownames(df2), df2[,i]), row.names = FALSE)
  datum <- datum[complete.cases(datum), ]
  species <- datum$V1
  char1 <- datum$V2
  names(char1) <- species
  
  for(j in c(1:ncol(df2))[-i]){
    # Was this combination of char1 and char2 already tested the other way around? Then skip
    if(nrow(subset(df.res, df.res$Char1 == colnames(df2)[j] & df.res$Char2 == colnames(df2)[i])) != 1){
      # reformat character 2  
      datum <- as.data.frame(cbind(rownames(df2), df2[,j]), row.names = FALSE)
      datum <- datum[complete.cases(datum), ]
      species <- datum$V1
      char2 <- datum$V2
      names(char2) <- species
      
      #pagels 1994 test of correlated evolution
      Pagel <- fitPagel(tree3, x=char1, y=char2, model = "ARD") # ARD model, also need to test ER
      # Pagel
      
      # Save p-value and log-likelyhood
      df.res.t <- data.frame('Char1' = colnames(df2)[i], 'Char2' = colnames(df2)[j], 'Pvalue' = Pagel[["P"]][1], 'LogLikelyhood' = Pagel[["lik.ratio"]][1])
      print(df.res.t)
      
      # Store in data.frame with all previous results
      if(exists('df.res') && is.data.frame(get('df.res'))){
        df.res <- rbind(df.res, df.res.t)
      }else{
        df.res <- df.res.t
      }
    }
  }
}

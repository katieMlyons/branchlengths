# setwd("/Users/miranda/Desktop/OpenTree/branchlengths")

library(geiger)

# TO RUN:  after loading all of the data and the function exchangeTaxa(), try ...
exchangeTaxa(tree, trait1, trait2) # runs on simulated data
exchangeTaxa(whales, with.br, randomwhales)

# library(phytools)

# Current code basically just groups taxa relative to other taxa of interest
# So, what we REALLY want is to take in a "big" tree (i.e., the synthetic
# tree or a subset of the synthetic tree), and a data.frame with comparative
# data and with branch length information.  That way, we can find sets of
# species that have both a species with comp data and with br data (if any
# exists).

#####################################################################
# Simulated example
tree <- sim.bdtree(b=1, d=0, n=50)
trait1 <- sample(tree$tip.label, 15)
trait2 <- sample(tree$tip.label, 5)

######################################################################
# Using the Cetacean data
whales <- read.tree("./data/pg_1927.tre") # pretending that it doesn't have branch lengths
bigtree <- read.tree("./data/phy_1428.phy") # they DO have branch lengths

# "trait1"
with.br <- whales$tip.label[which(whales$tip.label %in% bigtree$tip.label)]

# "trait2"
randomwhales <- sample(whales$tip.label, 10)


exchangeTaxa <- function(tree, sp.with.trait1, sp.with.trait2) {
  trait1 <- sp.with.trait1 # typically, comparative data
  trait2 <- sp.with.trait2 # typically, branch length data

  df <- data.frame(tree$tip.label, trait1=rep(0, length(tree$tip.label)), trait2=rep(0, length(tree$tip.label)), group=seq(from=1, to=length(tree$tip.label), by=1))
  df[which(df$tree.tip.label %in% trait1), "trait1"] <- 1 # sets each int.tax = 1 if that species is in the list of interesting taxa  
  df[which(df$tree.tip.label %in% trait2), "trait2"] <- 1 # sets each int.tax = 1 if that species is in the list of interesting taxa  
  
  for(i in trait1) {
    # get the immediate ancestral node for that interesting taxon
    # extract a subtree and get the list of species in that subtree
    anc <- tree$edge[which(tree$edge[,2]==(which(tree$tip.label==i))),1]
    sp.subset <- extract.clade(tree, anc)$tip.label
    
    # for that list of species, if there's one interesting taxon,
    # set all of the other groups = to the group of that interesting taxon
    if(sum(trait1 %in% sp.subset) == 1) {
      curr <- trait1[trait1 %in% sp.subset]
      group <- df[which(curr == df$tree.tip.label), "group"] # get the group for the interesting taxon
      df[which(df$tree.tip.label %in% sp.subset), "group"] <- group
    }
    
    # if there are more than one interesting taxa in there, stop
    if(sum(trait1 %in% sp.subset) > 1) {stop}
    
  }
  plot(tree, tip.color=df$group)
  
  # next up:  returning not just the df, but rather a list of taxa that
  # are exchangeable with trait2 (not trait1, because everything is
  # grouped by trait1)
  exchange.groups <- list()
  
  for(j in unique(df$group)) {
    # select the current group and extract the species in that group
    g <- df[which(df$tree.tip.label %in% j), "group"]
    sp <- as.character(df[which(df$group == j), "tree.tip.label"])
    
    # add that group to the list
    exchange.groups <- append(exchange.groups, list(sp))
  }
  return(exchange.groups)
#   return(list(exchange.groups, df))
}


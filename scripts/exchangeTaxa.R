library(geiger)
library(phytools)

# Current code basically just groups taxa relative to other taxa of interest
# So, what we REALLY want is to take in a "big" tree (i.e., the synthetic
# tree or a subset of the synthetic tree), and a data.frame with comparative
# data and with branch length information.  That way, we can find sets of
# species that have both a species with comp data and with br data (if any
# exists).


tree <- sim.bdtree(b=1, d=0, n=20)
int.tax <- sample(tree$tip.label, 10)
node <- fastMRCA(tree, int.tax[1], int.tax[2])

# identify the interesting taxa in the df
df <- data.frame(tree$tip.label, int.tax=rep(0, length(tree$tip.label)), group=seq(from=1, to=length(tree$tip.label), by=1))
df[which(df$tree.tip.label %in% int.tax), "int.tax"] <- 1 # sets each int.tax = 1 if that species is in the list of interesting taxa


# for each interesting taxon, go back one ancestor to find out if 
# there is another interesting taxon; if so, then stop
# otherwise, DROP THAT TIP, set its group equal to the group of the int.tax,
# and continue until it hits another interesting taxon
# must copy the tree first

group.tree <- tree

# very messy and ugly, but selects the species names that are not in the
# list of interesting taxa
not.int.tax <- as.character(df[which(!(df$tree.tip.label %in% int.tax)), "tree.tip.label"])

for(i in not.int.tax) {
  # get the immediate ancestral node for that interesting taxon
  # extract a subtree and get the list of species in that subtree
  anc <- tree$edge[which(tree$edge[,2]==(which(tree$tip.label==i))),1]
#   anc <- fastMRCA(group.tree, i, i)
#   which(tree$tip.label=='s3')
  
  sp.subset <- extract.clade(group.tree, anc)$tip.label
  
  # for that list of species, if there's one interesting taxon,
  # set all of the other groups = to the group of that interesting taxon
  if(sum(int.tax %in% sp.subset) == 1) {
    curr <- int.tax[int.tax %in% sp.subset]
    group <- df[which(curr == df$tree.tip.label), "group"] # get the group for the interesting taxon
    df[which(df$tree.tip.label %in% sp.subset), "group"] <- group
  }
  
  # if there are NO interesting taxa in there, set their group numbers 
  # equal
#   if(sum(int.tax %in% sp.subset) < 1) {
#     print(sp.subset)
#   }
  
  # if there are more than one interesting taxa in there, stop
  if(sum(int.tax %in% sp.subset) > 1) {
    stop
  }
}

plot(tree, tip.color=df$group)
# plot(tree, show.tip.label=FALSE)
# tiplabels(col=df$group, bg=df$int.tax)
df


# next two steps:
#   1.  fix the fastMRCA thing to get the ancestral node
#   2.  add in code for where two species are sister but don't have any
#       interesting taxa (? if it's necessary ?)



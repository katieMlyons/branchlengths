library(geiger)


# New approach!
# Instead of worrying about branch lengths, just identify a subset
# of the tree that you're interested in and collapse all other taxa
# into the same (appropriate) group.

### NOTE:  doesn't quite work correctly -- only groups together taxa that
# are SISTER to one of the "interesting taxa"
#  next step:  fix this problem!


tree <- sim.bdtree(b=1, d=0, n=10)
int.tax <- sample(tree$tip.label, 3)
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

for(i in int.tax) {
  print(i)
  
  # get the immediate ancestral node for that interesting taxon
  # extract a subtree and get the list of species in that subtree
  anc <- fastMRCA(group.tree, i, i)
  sp.subset <- extract.clade(group.tree, anc)$tip.label
  
  # for that list of species, if there's one interesting taxon,
  # set all of the other groups = to the group of that interesting taxon
  if(sum(int.tax %in% sp.subset) == 1) {
    curr <- int.tax[int.tax %in% sp.subset]
    group <- df[which(curr == df$tree.tip.label), "group"] # get the group for the interesting taxon
    df[which(df$tree.tip.label %in% sp.subset), "group"] <- group
  }
  
  # if there are more than one interesting taxa in there, stop
  if(sum(int.tax %in% sp.subset) != 1) {stop}  
  
}

plot(tree)
df

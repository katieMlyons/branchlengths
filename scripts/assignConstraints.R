# This set of functions generates constraints (calibrations) for a given tree, based on fossil occurrence data downloaded from the PBDB.
# For each node, the constraints are based on the first appearance of descendent nodes, rather than the first appearance of a given node. 
# This is to avoid assigning (potential) stem group taxa to the crown group. 
# In this respect, any constraints generated using this strategy will be conservative minimum estimates of divergence times.
# Presently, these functions only search for occurrence data based on genera that are present in the tree. 
# (In future implementations we may use higher taxonomy, which may increase the amount of occurrence data that could be used to inform calibrations.)
# Note that data downloaded from the PBDB should be treated with caution. (?Provide Canidae example?)


require(ape)
require(geiger)
require(phytools)
require(paleobioDB)

# 1. get the data (tree, occurrences)
# names="cetacea"
# pbdb.data<-pbdb_occurrences(limit="all",base_name=names,vocab="pbdb",show=c("phylo", "time", "ident"))

# tree <- read.tree("./data/pg_1927.phy")
# tree$tip.label <- gsub("\'", "", tree$tip.label)

# whale.tree <- tree
# tree <- extract.clade(tree, node=109)

# 2. get genus + age related info
# g <- getGenera(tree,pbdb.data)

# 3. assign calibrations to each node
# calibrations<-getCalibrations(g)

# 4. removed nested calibrations
# fixed.constraints<-filter.constraints(calibrations,tree)

# 5. assignning daughter nodes to calibrated parent nodes
# calibrated.nodes<-fill.daughters(fixed.constraints,tree)

# Obtain age related information for the descendents of each node.
getGenera <- function(tree, pbdb.data, upper.bound=100) {
# This function takes a tree + pbdb occurrence data.
# For each (active) node it gets the descendent genera, extracts the occ. data for each genus
# works out the oldest age of each genus, and assign this age to the active node.
# The function outputs a dataframe containing this information.
# Note this function will only work for trees containing tips in the format Genus_species
  
  # eliminare apostrophes
  tree$tip.label <- gsub("\'", "", tree$tip.label)
  
  # for each node this block finds its two desecendent daughter nodes 
  # (and all their descendents)
  genera<-data.frame(node=numeric(), genus=character(), age=numeric())
  for (i in (length(tree$tip.label)+1):(length(tree$tip.label)*2-1) ) {
      row=which(tree$edge[,1]==i)
      descendents=tree$edge[,2][row]
      dec1 <- getDescendants(tree, descendents[1])
      dec1 <- dec1[which(dec1 <= length(tree$tip.label))]
      dec2 <- getDescendants(tree, descendents[2])
      dec2 <- dec2[which(dec2 <= length(tree$tip.label))]
      
      # input the node number and all descendent genera into a data frame
      s=strsplit(tree$tip.label[c(dec1, dec2)], "_")
      for(j in 1:length(s)) {
        genera<-rbind(genera,data.frame(node=i,genus=s[[j]][1],age=NA))
      }
    }
  
  # the following obtains the youngest secure age of the oldest fossil for each genus using the pbdb occurrence data
  for(k in as.character(unique(genera$genus))) {
    age <- pbdb.genus.min(pbdb.data, k,upper.bound)
    genera[which(k == genera$genus), "age"] <- age
  }
  
  return(unique(genera))
} # end of function
  
# g <- getGenera(tree,pbdb.data)

# ------ 
# Identify calibrations
getCalibrations<-function(genera,upper.bound) {
  # this function generates a calibration matrix (node and corresponding age constraints)
  # based on pbdb occurrence data
  
  # assign that date to the ancestral node
  calibrations <- data.frame(node=unique(g$node), age=NA)

  for(n in 1:nrow(calibrations)) {
    # pick out the current active node
    active.node <- calibrations[n, "node"]
  
    # find the number of times it occurs in the "g" data.frame (if it occurs more than once,
    # it had more than one descendent genus, meaning that we can calibrate it)
    row=which(g$node==active.node)
  
    # when the node has more than one descendent genus:
    if( length(row) > 1) {
      # get the list of ages for those genera
      ages <- g[which(g$node == active.node), ]
    
      # as long as there's at least one age, find the max and keep
      if(nrow(ages) > 0) {
#         print(ages)
        max.age <- max(na.omit(ages$age))
        calibrations[n, "age"] <- max.age # stores max age in calibrations
      } 
    }
  }
return(calibrations)
# EOF  
}

# calibrations<-getCalibrations(g)

# ------ 

# Eliminate nested calibrations
filter.constraints<-function(calibrations,tree) {
  # the purpose of this function is to eliminate duplicate (nested) calibrations
  # e.g. ancestor descendent pairs that share the same calibration (because the first appearance in the fossil record is the same for both nodes)
  # The function returns the node - calibration pair for the youngest node only.

fixed.constraints<-data.frame(node=numeric(), age=numeric())

for(c in 1:nrow(calibrations)) {
  
  # select the current node
  active.node <- calibrations[c, "node"]
  active.minima<- calibrations[c, "age"]
  
  # check if calibration data is numerical
  if(is.finite(active.minima)) {
    
    # find all the descendants of that node
    decs<-getDescendants(tree, active.node)
    
    #decs.min=list()
    decs.min=0
    
    # for each descendent:
    for(d in 1:length(decs)) {
      
      # check if they occur in the calibration table
      if(decs[d] %in% calibrations$node) {
        row=which(calibrations$node==decs[d])
        e=calibrations$age[row]
        
        # check if calibration data is numerical
        if(is.finite(e)) {
          
          # add it to the list of constraints available for descendent nodes
          decs.min=c(decs.min,e)
          
        }      
      }
    }
    
    # if descendent calibrations are available
    if(length(decs.min) > 1) {
      
      # calculate the age of the oldest
      decs.min.max=max(decs.min)
      
      # only apply the calibration to the current (active) node 
      # if the max descendent calibration is NOT > = the active calibration (active.minima)
      if (!decs.min.max >= active.minima){
        fixed.constraints<-rbind(fixed.constraints,data.frame(node=active.node,age=active.minima))        
      }
    }
    # if no descendent calibrations are available
    # apply the current calibration (active.minima) to the current (active) node
    else {
      fixed.constraints<-rbind(fixed.constraints,data.frame(node=active.node,age=active.minima))
    }
  }  
}

return(fixed.constraints)

# EOF

}

# fixed.constraints<-filter.constraints(calibrations,tree)

# ------ 
# Create the table required for congruify
fill.daughters<-function(fixed.constraints,tree) {
# this function simply creates a table containing (1) node ID, (2) calibration age, (3) descendent nodes IDs, and 
# (4) respentative tip labels representing the divergence. Note this table will contain non-nested calibations.

# add the descendent node IDs to the fixed.constraints table
fixed.constraints<-cbind(fixed.constraints,child1=0,child2=0,tip1=0,tip2=0) 

for(i in 1:nrow(fixed.constraints)) {
  active.node <- fixed.constraints[i, "node"]
  row=which(tree$edge[,1]==active.node)
  descendents=tree$edge[,2][row]
  # get descendent node/edge IDs
  fixed.constraints$child1[i]=descendents[1]
  fixed.constraints$child2[i]=descendents[2]
  
  # get representative tip IDs & lables
  decs<-getDescendants(tree, descendents[1])
  
  if(length(decs) > 1) {
    rand.tip.no=(sample(decs,1))    
  }
  else {
    rand.tip.no=decs[1]
  }
  rand.tip1=tree$tip.label[rand.tip.no]
  
  decs<-getDescendants(tree, descendents[2])
  if(length(decs) > 1) {
    rand.tip.no=(sample(decs,1))    
  }
  else {
    rand.tip.no=decs[1]
  }
  rand.tip2=tree$tip.label[rand.tip.no]
  
  fixed.constraints$tip1[i]=rand.tip1
  fixed.constraints$tip2[i]=rand.tip2
  
}
return(fixed.constraints)
# EOF
}

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

getGenera <- function(tree, pbdb.data) {
# this function takes a tree + pbdb occurrence data
# for each (active) node it gets the descendent genera, extracts the occ. data for each genus
# works out the oldest age of each genus, and reassign to the active node
# the function will only work for trees containing tips of the format 
# Genus_species  
  
  # eliminare apostrophes
  tree$tip.label <- gsub("\'", "", tree$tip.label)
  
  # for each node this block finds its desecendents
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
  
  # this obtains the youngest secure age of the oldest fossil for each genus using the pbdb occurrence data
  for(k in as.character(unique(genera$genus))) {
    age <- pbdb.genus.min(pbdb.data, k)
    genera[which(k == genera$genus), "age"] <- age
  }
  
  return(unique(genera))
} # end of function
  
# g <- getGenera(tree,pbdb.data)

####

# this function will generate a calibration matrix
# and populate it using constraints based on pbdb occurrence data
getCalibrations<-function(genera,upper.bound) {

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

# filter.constraints
# the purpose of this function is to eliminate duplicate (nested) calibrations

filter.constraints<-function(calibrations,tree) {

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
# the function fill.daughters
# inserts daughter node IDs into the (non-nested) calbration table

fill.daughters<-function(fixed.constraints,tree) {

# add the descendent node IDs to the fixed.constraints table
fixed.constraints<-cbind(fixed.constraints,child1=0,child2=0) 

for(i in 1:nrow(fixed.constraints)) {
  active.node <- fixed.constraints[i, "node"]
  row=which(tree$edge[,1]==active.node)
  descendents=tree$edge[,2][row]
  fixed.constraints$child1[i]=descendents[1]
  fixed.constraints$child2[i]=descendents[2]  
}
return(fixed.constraints)
# EOF
}

## Functions
content <- httr:::content # does this mean this is necessary? require(httr)

## Get all studies with these ottids
findsourcetrees <- function(otts){
  ##start with a list of ottids
  findstudies <- function(otts){
    res <- rotl::studies_find_trees(property="ot:ottId", value=otts,exact = TRUE)
    res <- fromJSON(content(res, "text"))$matched_studies
  }
  gettrees <- function(x, studies){
    matched_trees <- studies[[x]]
    treeId <- sapply(matched_trees, function(x) x$matched_trees[[1]]$nexson_id)
    studyId <- sapply(matched_trees, function(x) x['ot:studyId'][[1]])
    df <- cbind(studyId, treeId)
  }
  
  studies <- lapply(otts, function(x) findstudies(x))
  table <- do.call(rbind, lapply(1:length(studies), gettrees, studies=studies))
  counts <- table(paste(table[,1], table[,2], sep="_"))
  table <- unique(table)
  counts <- counts[paste(table[,1], table[,2], sep="_")]
  list(counts=counts, table = table)
}

## Download the trees of a specific study and tree
.getStudyTrees <- function(study, phy, format="newick"){
  res <- httr:::GET(paste(rotl:::otl_url(),"/", rotl:::otl_version(), "/study/", study, "/tree/", phy, sep=""))
  stringRes <- content(res,"text")
  jsRes <- fromJSON(stringRes)
  hasBranches <- length(grep("@length", stringRes))>0
  branchLengthDescription <- jsRes[[1]]['^ot:branchLengthDescription']
  branchLengthMode <- jsRes[[1]]['^ot:branchLengthMode']
  branchLengthTimeUnit <- jsRes[[1]]['^ot:branchLengthTimeUnit']
  if(hasBranches){
    nwTree <- httr:::GET(paste(rotl:::otl_url(),"/", rotl:::otl_version(), "/study/", study, "/tree/", phy,".tre", sep=""))
    tree_string <- httr:::content(nwTree, "text")
    tree_string <- gsub( "\\[.*?\\]", "", tree_string)
    tree_string <- gsub(" ", "_", tree_string)
    tree <- read.tree(text=tree_string)
  } else {
    return(NULL)
  }
  rm(jsRes)
  return(list=c(tree=tree, ntips= length(tree$tip.label), branchLengthDescription=branchLengthDescription,
                branchLengthMode=branchLengthMode, branchLengthTimeUnit=branchLengthTimeUnit))
}
getStudyTrees <- function(studyTrees){  
  res <- lapply(which(studyTrees$counts>2), function(x) .getStudyTrees(studyTrees$table[x,1], studyTrees$table[x,2], format="newick"))
  res <- res[which(!(sapply(res, is.null)))]
  summary <- data.frame(ntips = sapply(res, function(x) x$ntips),
                        matching = studyTrees$counts[names(res)],
                        timeUnits = sapply(res, function(x) x$branchLengthTimeUnit[[1]]),
                        description = sapply(res, function(x) x$branchLengthDescription[[1]]),
                        mode = sapply(res, function(x) x$branchLengthMode[[1]]))
  return(list(result=res, summary=summary))
}

extractTree <- function(i, possibleTrees){
  res <- possibleTrees$res
  tree <- res[[i]][[1]]
  tree$tip.label <- gsub("'", "", tree$tip.label, fixed=TRUE) ## tidy up this tree, by removing apostrophes
  return(tree)
}

### Search database for calibration points
## Find all trees with calibration times
sourceTreeCalibrations <- function(tree, taxalist){
  newTaxa <- tree$tip.label[!(tree$tip.label %in% taxalist)]
  newTaxa <- gsub("_", " ", newTaxa)
  newTaxa <- gsub(" sp[.]", "", newTaxa)
  newTaxa <- gsub(" spp[.]", "", newTaxa)

  ## faster to turn off approximate matching BUT may return NAs
  newTax <- rotl::tnrs_match_names(newTaxa, do_approximate_matching = FALSE)
  mmTaxa <- which(is.na(newTax$ott_id))
  if(length(mmTaxa)>0){
    mmTax <- rotl::tnrs_match_names(newTaxa[mmTaxa], do_approximate_matching = TRUE)
    newTax[mmTaxa,] <- mmTax
  }
  ## Find all timetrees
  tttaxa <- newTax$ott_id
  ttsource <- findsourcetrees(tttaxa)
  ## must have at least 2 matching tips
  ttres <- lapply(which(ttsource$counts>2), function(x) .getStudyTrees(ttsource$table[x,1], ttsource$table[x,2], format="newick"))
  ttres <- ttres[which(!(sapply(ttres, is.null)))]
  ## must contain branch lengths
  timetrees <- ttres[which(sapply(ttres ,function(x) x$branchLengthMode == "ot:time"))]
  ttsummary <- data.frame(ntips = sapply(timetrees, function(x) x$ntips),
                        matching = ttsource$counts[names(timetrees)],
                        timeUnits = sapply(timetrees, function(x) x$branchLengthTimeUnit[[1]]),
                        description = sapply(timetrees, function(x) x$branchLengthDescription[[1]]),
                        mode = sapply(timetrees, function(x) x$branchLengthMode[[1]])
  )
  ## get all timetrees, drop non-matching tips
  tmp <- lapply(timetrees, function(x) x$tree)                  
  dropAbsent <- function(x, taxalist){
      tr <- x$tree
      tr$tip.label <- gsub("\'", "", tr$tip.label)
    tipL <- tr$tip.label
    toDrop <- tipL[!(tipL %in% taxalist)]
    tr <- drop.tip(tr, toDrop)
    return(tr)
  }
  #tmp[[1]]$tip.label
  dropTipTrees <- lapply(timetrees, function(x) dropAbsent(x, gsub(" ", "_", newTaxa)))
  ## in case any edge lengths are NA, remove those trees
  dropna <- function(tlist){
    nalist <- lapply(dropTipTrees, function(x) any(is.na(x$edge.length)))
    trues <- grep(TRUE, nalist)
    tlist[trues] <- NULL
    tlist
  }
  cleaned <- dropna(dropTipTrees)
  return(list(summary=ttsummary, cleanTrees=cleaned))
}


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
}

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
  fixed.constraints<-cbind(fixed.constraints,child1=0,child2=0,tip1=0,tip2=0) 
  
  for(i in 1:nrow(fixed.constraints)) {
    active.node <- fixed.constraints[i, "node"]
    row=which(tree$edge[,1]==active.node)
    descendents=tree$edge[,2][row]
    # get descendent node/edge IDs
    fixed.constraints$child1[i]=descendents[1]
    fixed.constraints$child2[i]=descendents[2]
    
    # get representative tip IDs
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


pbdb.genus.min<-function(pbdb.data,genera, upper.bound=100)
  # note the upper bound is fixed but you may with to change this depending on the group
  {
  
  if(!length(pbdb.data)==0) { # if the data frame is not empty
    
    # extract columns for genus only
    rows=which(pbdb.data$genus_name==genera)    
    late.ages=pbdb.data$late_age[rows] # colunm contains: late bound of the geologic time range associated with this occurrence (in Ma)
    late.ages=sort(late.ages) # sort a -> b
    
    ### the following is a quick solution
    ### for eliminating dates that are implausibly old
    count=0
    
    while(count==0) {
      
      if (length(late.ages) > 0) { # if there are available dates
        
        max.age=max(late.ages) # define the greatest
        
        if(max.age >=upper.bound) { # if this is greater than the specified upper limit
          late.ages<-late.ages[-length(late.ages)] # remove the last element of this list
        }
        else {
          
          count=1  
        }
      }
      # this would mean that there were no suitable dates available in the list
      else {
        max.age=0
        count=1
      }	
    }
    ###
  }
  
  else {
    max.age=0;
  }
  
  if(max.age!=0) {
    return(max.age) # this represents the minimum maximum age!
  }
  else {
    return(NA)
    #     print("I'm sorry, no calibration info available!")
    
  }
  #EOF
}

getSourceTreeCalibrations <- function(cal){
  tmp <- lapply(cal$cleanTrees, function(x) congruify.phylo(x, target, scale=NA))
  calibrations <- lapply(tmp, function(x) x$calibrations)
  calibrations <- do.call(rbind, calibrations)
  return(calibrations)
}

replaceHashes <- function(target, table){
  scion <- target
  taxa <- unique(unlist(table[,4:5]))
  scion <- hashes.phylo(scion, taxa)
  hashes <- lapply(1:nrow(table), function(x) scion$hash[getMRCA(target, c(table[x,4], table[x,5]))])
  table[,1] <- unlist(hashes)
  return(list(scion=scion, table=table))
}
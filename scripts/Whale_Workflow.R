## # Complete whale analysis workflow:

#################################

## ## Part 1. Set the working directory, load required libraries and functions
# setwd("~/repos/OThackathon/branchlengths/scripts/")
# setwd("/Users/miranda/Desktop/OpenTree/branchlengths/scripts/")

require(devtools)
install_github("rotl", user="fmichonneau")
# require(rotl)
# require(rjson)
# require(geiger)
# require(httr)

#################################

## ## Part 2. Load required input files
## The file CetaceanCompData.csv contains a list of taxa for which we have compartive data.
table <- read.csv("../data/CetaceanCompData.csv")
bldatabase <- read.csv("../data/branch_length_studies.csv", header=FALSE) # where will the users obtain this?

content <- httr:::content # does this mean this is necessary? require(httr)

## TNRS
#require(devtools) # delete
#install_github("rotl", user="fmichonneau") # delete
require(rotl) # delete
require(rjson) # delete
require(geiger) # delete

#################################

## ## Part 3. Download source trees from the Open Tree of Life database.

# Step 1. Generate a list of taxa for which we have comparive data.
taxalist <- table$genspec

# Step 2. what are the following two lines doing?
tax <- rotl::tnrs_match_names(as.character(taxalist))

#tax <- rotl::tnrs_match_names(taxalist)

tax <- tnrs_match_names(as.character(taxalist))

# Step 3. Gets ottIDs for the comparative data.
ottids <- as.character(tax$ott_id) 

# Step 4. Download all studies (can we call these source trees?) that contain these ottIDs.
# (the function findsourcetrees has been moved to the bottom of this file)
studyTrees <- findsourcetrees(ottids)      

## Step 5. Filter the (source?) trees to obtain a set that have branch lengths.
# (the function getStudyTrees has been moved to the bottom of this file)

#tmp <- fromJSON(content(rotl:::get_study(study = "pg_1816"),"text")) # delete
#tmp[[6]] # delete

# the following block looks like it could become a seperate function

res <- lapply(which(studyTrees$counts>2), function(x) getStudyTrees(studyTrees$table[x,1], studyTrees$table[x,2], format="newick"))
res <- res[which(!(sapply(res, is.null)))]
summary <- data.frame(ntips = sapply(res, function(x) x$ntips),
                      matching = studyTrees$counts[names(res)],
                      timeUnits = sapply(res, function(x) x$branchLengthTimeUnit[[1]]),
                      description = sapply(res, function(x) x$branchLengthDescription[[1]]),
                      mode = sapply(res, function(x) x$branchLengthMode[[1]])
)

# Step 6. Get tree with the largest number of taxa.
##Chosen tree and study: 1816 and tree3671
whichTree <- 2
tree <- res[[whichTree]][[1]]
tree$tip.label <- gsub("'", "", tree$tip.label, fixed=TRUE) ## tidy up this tree, by removing apostrophes

#################################

## ## Part 4. Find all time calibrated trees.

## # Find all trees with calibration times
#require(foreach)
#require(doParallel)
#registerDoParallel(cores=ncores)
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
ttres <- lapply(which(ttsource$counts>2), function(x) getStudyTrees(ttsource$table[x,1], ttsource$table[x,2], format="newick"))
ttres <- ttres[which(!(sapply(ttres, is.null)))]
## must contain branch lengths
timetrees <- ttres[which(sapply(ttres ,function(x) x$branchLengthMode == "ot:time"))]
ttsummary <- data.frame(ntips = sapply(timetrees, function(x) x$ntips),
                        matching = ttsource$counts[names(timetrees)],
                        timeUnits = sapply(timetrees, function(x) x$branchLengthTimeUnit[[1]]),
                        description = sapply(timetrees, function(x) x$branchLengthDescription[[1]]),
                        mode = sapply(timetrees, function(x) x$branchLengthMode[[1]])
)
ttsummary

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
tmp[[1]]$tip.label

dropTipTrees <- lapply(timetrees, function(x) dropAbsent(x, gsub(" ", "_", newTaxa)))
## in case any edge lengths are NA, remove those trees
dropna <- function(tlist){
  nalist <- lapply(dropTipTrees, function(x) any(is.na(x$edge.length)))
  trues <- grep(TRUE, nalist)
  tlist[trues] <- NULL
  tlist
}
cleaned <- dropna(dropTipTrees)
#################################
## cleaned is the list of timetrees overlapping with the chosen source tree!

## # Get taxonomy
tree$tip.label <- gsub("\'", "", tree$tip.label)
taxonomy  <- gbresolve(tree)
plot(taxonomy$phy, show.node=TRUE, type="f", cex=0.5)

#################################

## ## Part 5. Download calibration data from the Paleobiology Database.
# This section describes how to utilize fossil information available on the Paleobiology Database (http://paleobiodb.org/) via the PBDB API.

# The following requires the PaleobioDB R library (https://github.com/paleobiodb/paleobiodb_utilities).
require(paleobioDB)

# Step 1. Specify the most inclusive group(s) of taxa represent in the tree. This may also be a list of taxa.
names="cetacea"
# Step 2. Download the data from PBDB. More information about options from downloading data can be obtained here (http://cran.r-project.org/web/packages/paleobioDB/paleobioDB.pdf) and here (http://paleobiodb.org/data1.1/occs/list).
pbdb.data<-pbdb_occurrences(limit="all", base_name=names, vocab="pbdb", show=c("phylo", "time", "ident"))

# Step 3. Get age related information from the data downloaded from PBDB for each genus in the tree.
# The function getGenera finds the oldest fossil occurrence for each genera and assigns it the youngest secure age of the associated fossil occurrence.
g <- getGenera(tree, pbdb.data, upper.bound=100)

# Step 4. Assign calibrations to each node.
# The function getCalibrations assigns calibration information to as many nodes in the tree as possible.
calibrations<-getCalibrations(g) # produces a lot of -Inf for calibrations and consequently warnings; however, still runs

# Step 4. Eliminate nested (redundant) calibrations. This is neccesary because we're using fixed calibrations.
fixed.constraints<-filter.constraints(calibrations,tree)

# Step 5. Identify the daughters (and representative tip labels) assiciated with each calibrated parent node.
calibrated.nodes<-fill.daughters(fixed.constraints,tree)

#################################

## ## Part 7. Build final timetree.

## congruify
res <- congruify.phylo(reference, target)





#################################

## ## Part 6. Identify a set of exchangable taxa.

# At this stage we have a set of species (on a tree) with comparative data,
# and a set of species (on a different tree) with branch lengths, and a tree
# (the scaffold/synthetic OpenTree) that includes most/all of the species that are in the 
# first two trees.  We want to use the synthetic tree to match up names on
# comparative tree with names on the time tree to figure out if any of the
# comparative tree species can be exchanged with species on the time tree.

# Step 1. Read the tree containing the species for which we have comparative data.
# comparative.tree <- read.tree("../data/Induced_whaletree.phy")
time.tree <- tree # rename the tree to time.tree to clarify its purpose here

#  Step 2:  get the synthetic tree by getting the subtree for the root of the comparative tree
mrca <- tol_mrca(ott_ids=ottids)$mrca_node_id
scaffold <- tol_subtree(node_id=mrca)
scaffold$tip.label <- gsub("_ott\\d+", "", scaffold$tip.label) # remove ottID after the tip labels
# comparative.tree$tip.label <- gsub("_ott\\d+", "", comparative.tree$tip.label) # remove ottID after the tip labels
# time.tree <- tree

# Step 3. Identify a set of exchangable taxa.
# (the function exchangeTaxa has been moved to the bottom of this file)
taxa.df <- findExchangeableGroups(scaffold, taxalist)
exchange.groups <- listExchangeableTaxa(taxa.df)
tree <- exchangeTaxa(time.tree, taxalist, exchange.groups)  


#################################





#################################
# Functions that should be in a seperate file

# Function findsourcetrees
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

# Function getStudyTrees
getStudyTrees <- function(study, phy, format="newick"){
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


##########################################################
###  Getting taxa you care about onto a tree you have time information for
###
###  Sometimes you have comparative data for species that you're interested in,
###  and access to a time tree that has species you're not as excited about. Using
###  a larger, more inclusive scaffold tree (such as a portion of the synthetic
###  tree derived from OpenTree), you can match species for which you have
###  comparative data to species for which you have time data, based on their
###  phylogenetic relationships: if a species with comparative data is closely related
###  to one with time data, you can exchange those taxa to create a version of the
###  time tree that includes more of the species you're interested in.
##########################################################
# Function getAnc
getAnc <- function(tree, node) {
  if(is.character(node)) {
    node <- which(tree$tip.label == node)
  }
  
  row <- which(tree$edge[,2]==node)
  ancestor <- tree$edge[,1][row]
  return(ancestor)
}

# Function findExchangeableGroups
findExchangeableGroups <- function(tree, species) {
  trait <- species # typically, comparative data
  
  # setting up the data frame that contains the list of species, whether they
  # have comparative data, and whether they have branch length data
  df <- data.frame(tree$tip.label, "trait"=rep(0, length(tree$tip.label)), "group"=seq(from=1, to=length(tree$tip.label), by=1))
  df[which(df$tree.tip.label %in% trait), "trait"] <- 1 # sets each int.tax = 1 if that species is in the list of interesting taxa
  
  # list of comparative species and their exchangeable taxa
  exchange.groups <- list()
  
  # loop through each species with comparative trait data (in trait1)
  for(i in trait) {
    
    # if and only if the species is present in the big tree:
    if(i %in% tree$tip.label) {
      
      count <- 1 # count of the number of comparative species in the current species set
      node <- i # node to get the ancestor for; can be a node number or species name
      
      while(count < 2) {
        anc <- getAnc(tree, node) # get the ancestral node
        sp <- extract.clade(tree, anc)$tip.label # extract the subtree and get sp names
        
        # update the count and the node/tip
        count <- sum(sp %in% trait)
        node <- anc
        
        if(count < 2) {
          curr <- sp[sp %in% trait]
          group <- df[which(curr == df$tree.tip.label), "group"] # get the group for the interesting taxon
          df[which(df$tree.tip.label %in% sp), "group"] <- group
        }
      } # end of while loop
    } # end of what to do if the species isn't in the big tree 
  } # end of looping through species with comparative data
  
  return(df)
  
}

# Function listExchangeableTaxa
listExchangeableTaxa <- function(df) {
  comp.sp <- df[which(df$trait == 1), ]
  exchange.groups <- list()
  for(j in 1:nrow(comp.sp)) {
    # get the group number for that species
    group.num <- comp.sp[j, "group"]
    
    # extract all the species that share that group number
    e.taxa <- as.character(df[which(df$group == group.num), "tree.tip.label"])
    
    # add those species to the exchange.groups object
    exchange.groups <- append(exchange.groups, list(e.taxa))
  }
  
  return(exchange.groups)
}

# Function exchangeTaxa
exchangeTaxa <- function(time.tree, comparative, exchange.groups) {
  time.sp <- time.tree$tip.label
  exchange.taxa.df <- data.frame("branch.length"=NA, "comparative.exchange"=NA)
  
  # for each list of exchangeable taxa
  for(group in exchange.groups) {
    # if and only if there's a comparative species AND a branch length species in the group
    if(sum(group %in% time.sp) > 0) {
      br.sp <- group[which(group %in% time.sp)] # should check if the number of species with branch length info is more than 1
      comparative.sp <- group[which(group %in% comparative)]
      
      exchange.taxa.df <- rbind(exchange.taxa.df, c(br.sp, comparative.sp))
    } # end of if-statement
  } # end of looping through exchange groups
  
  exchange.taxa.df <- exchange.taxa.df[2:nrow(exchange.taxa.df),] # removes the first dummy row
  
  
  # using the list of branch length species, subset the branch length tree
  for(s in time.tree$tip.label) {
    if(!(s %in% exchange.taxa.df$branch.length)) {
      time.tree <- drop.tip(time.tree, s)
    }
  }
  
  # switch out the names
  for(t in 1:length(time.tree$tip.label)) {
    tip <- time.tree$tip.label[t]
    row <- which(exchange.taxa.df$branch.length == tip)
    new.name <- exchange.taxa.df[row, "comparative.exchange"]
    time.tree$tip.label[t] <- new.name
  }
  
  return(time.tree)
}

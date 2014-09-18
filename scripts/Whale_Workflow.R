## Complete whale analysis workflow:
## Setup
# setwd("~/repos/OThackathon/branchlengths/scripts/")
setwd("/Users/miranda/Desktop/OpenTree/branchlengths/scripts/")
## inputs:
table <- read.csv("../data/CetaceanCompData.csv")
bldatabase <- read.csv("../data/branch_length_studies.csv", header=FALSE)
content <- httr:::content


## TNRS
#require(devtools)
#install_github("rotl", user="fmichonneau")
require(rotl)
require(rjson)
require(geiger)
taxalist <- table$genspec
tax <- rotl::tnrs_match_names(taxalist)
tax <- tnrs_match_names(as.character(taxalist))

ottids <- as.character(tax$ott_id) # gets ottID for the comparative data

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
studyTrees <- findsourcetrees(ottids)      

## Filter down to trees that have branch lengths
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

#tmp <- fromJSON(content(rotl:::get_study(study = "pg_1816"),"text"))
#tmp[[6]]


res <- lapply(which(studyTrees$counts>2), function(x) getStudyTrees(studyTrees$table[x,1], studyTrees$table[x,2], format="newick"))
res <- res[which(!(sapply(res, is.null)))]
summary <- data.frame(ntips = sapply(res, function(x) x$ntips),
                      matching = studyTrees$counts[names(res)],
                      timeUnits = sapply(res, function(x) x$branchLengthTimeUnit[[1]]),
                      description = sapply(res, function(x) x$branchLengthDescription[[1]]),
                      mode = sapply(res, function(x) x$branchLengthMode[[1]])
)

## Get tree with most taxa
##Chosen tree and study: 1816 and tree3671
whichTree <- 2
tree <- res[[whichTree]][[1]]
tree$tip.label <- gsub("'", "", tree$tip.label, fixed=TRUE)






#############################################################################
#############################################################################
#############################################################################
# Exchange taxa:  we have a set of species (on a tree) with comparative data,
# and a set of species (on a different tree) with branch lengths, and a tree
# (the synthetic tree) that includes most/all of the species that are in the 
# first two trees.  We want to use the synthetic tree to match up names on
# comparative tree with names on the time tree to figure out if any of the
# comparative tree species can be exchanged with species on the time tree.

comparative.tree <- read.tree("../data/Induced_whaletree.phy")
time.tree <- tree

###  Step 1:  get the synthetic tree by getting the subtree for the root of
# the comparative tree
mrca <- tol_mrca(ott_ids=ottids)$mrca_node_id
subtree <- tol_subtree(node_id=mrca)
subtree$tip.label <- gsub("_ott\\d+", "", subtree$tip.label) # remove ottID after the tip labels
comparative.tree$tip.label <- gsub("_ott\\d+", "", comparative.tree$tip.label) # remove ottID after the tip labels




exchangeTaxa <- function(tree, sp.with.trait1, sp.with.trait2) {
  trait1 <- sp.with.trait1 # typically, comparative data
  trait2 <- sp.with.trait2 # typically, branch length data

  # setting up the data frame that contains the list of species, whether they
  # have comparative data, and whether they have branch length data
  df <- data.frame(tree$tip.label, trait1=rep(0, length(tree$tip.label)), trait2=rep(0, length(tree$tip.label)), group=seq(from=1, to=length(tree$tip.label), by=1))
  df[which(df$tree.tip.label %in% trait1), "trait1"] <- 1 # sets each int.tax = 1 if that species is in the list of interesting taxa
  df[which(df$tree.tip.label %in% trait2), "trait2"] <- 1 # sets each int.tax = 1 if that species is in the list of interesting taxa

  # list of comparative species and their exchangeable taxa
  exchange.groups <- list()
  
  # loop through each species with comparative trait data (in trait1)
  for(i in trait1) {
    
    # if and only if the species is present in the big tree:
    if(i %in% tree$tip.label) {
      
      # get the immediate ancestral node for that interesting taxon
      # extract a subtree and get the list of species in that subtree
      anc <- tree$edge[which(tree$edge[,2]==(which(tree$tip.label==i))),1]
      sp.subset <- extract.clade(tree, anc)$tip.label
      
      # for that subset of species, if there's one interesting taxon,
      # set all of the other groups = to the group of that interesting taxon
      if(sum(trait1 %in% sp.subset) == 1) {
        curr <- trait1[trait1 %in% sp.subset]
        group <- df[which(curr == df$tree.tip.label), "group"] # get the group for the interesting taxon
        df[which(df$tree.tip.label %in% sp.subset), "group"] <- group
      }
      
      # if there are more than one interesting taxa in there, stop
      if(sum(trait1 %in% sp.subset) > 1) {
        stop
      }
      
    } # end of what to do if the species isn't in the big tree 
    
  } # end of looping through species with comparative data

  # next up: return the list of taxa that are exchangeable with each of
  # the comparative species
  comp.sp <- df[which(df$trait1 == 1), ]
  for(j in 1:nrow(comp.sp)) {
    # get the group number for that species
    group.num <- comp.sp[j, "group"]
    
    print(group.num)
    
    # extract all the species that share that group number
    e.taxa <- as.character(df[which(df$group == group.num), "tree.tip.label"])
    
    # add those species to the exchange.groups object
    exchange.groups <- append(exchange.groups, list(e.taxa))
  }

  return(exchange.groups)

}


sim.taxa <- exchangeTaxa(subtree, comparative.tree$tip.label, time.tree$tip.label)





#############################################################################
#############################################################################
#############################################################################



## Find all trees with calibration times
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

## Get taxonomy
tree$tip.label <- gsub("\'", "", tree$tip.label)
taxonomy  <- gbresolve(tree)
plot(taxonomy$phy, show.node=TRUE, type="f", cex=0.5)

## Get calibrations


## Build final timetree


## congruify
res <- congruify.phylo(reference, target)



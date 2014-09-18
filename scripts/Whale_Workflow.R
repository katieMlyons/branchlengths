## Complete whale analysis workflow:
## Setup
setwd("~/repos/OThackathon/branchlengths/scripts/")
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
ottids <- as.character(tax$ott_id)

## Get all studies with these ottids
findsourcetrees <- function(otts){
  ##start with a list of ottids
  findstudies <- function(otts){
    res <- rotl::studies_find_trees(property="ot:ottId", value=otts)
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

## Find all trees with calibration times
#require(foreach)
#require(doParallel)
#registerDoParallel(cores=ncores)
newTaxa <- tree$tip.label[!(tree$tip.label %in% taxalist)]
newTaxa <- gsub("_", " ", newTaxa)
newTax <- rotl::tnrs_match_names(newTaxa)






## Find all timetrees


## Replace taxa that can be replaced with exchangeable taxa 

## Get taxonomy
tree$tip.label <- gsub("\'", "", tree$tip.label)
taxonomy  <- gbresolve(tree)
plot(taxonomy$phy, show.node=TRUE, type="f", cex=0.5)

## Get calibrations


## Build final timetree


## congruify
res <- congruify.phylo(reference, target)



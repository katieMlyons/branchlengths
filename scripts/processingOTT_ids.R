## Get required packages and functions
require(rjson)
require(devtools)
install_github("rotl", user="fmichonneau")
content <- httr:::content

## Get ottids from the data table
dat <- read.csv("../data/CetaceanCompData.csv")
taxalist <- dat$genspec
tax <- rotl::tnrs_match_names(taxalist)
otts <- tax$ott_id

## Get the least inclusive ancestor
res <- rotl:::taxonomy_lica(ott_ids = otts[1:2])
res <- fromJSON(content(res, "text")) # parsing JSON
cladename <- res$lica$'ot:ottTaxonName'
cladename

## Run this to map the taxonomy onto the tree:
require(geiger)
require(ncbit)
tree <- read.tree("../data/pg_1927.phy")
tree$tip.label <- gsub("\'", "", tree$tip.label)
taxonomy  <- gbresolve(tree) # find the complete taxonomy from ncbi
plot(taxonomy$phy, show.node=TRUE, type="f", cex=0.5)

## See the taxonomy here:
taxonomy$tax

## You'll need to match these to the ottids to find the lica
require(phytools)
getDescTips <- function(tree){
  N <- length(tree$tip.label)
  internal_nodes <- (N+1):(N+tree$Nnode)
  descendants <- lapply(internal_nodes, function(x) getDescendants(tree, x))
  descendant_tips <- lapply(descendants, function(x) x[x < N])
  descendant_tips <- lapply(descendant_tips, function(x) tree$tip.label[x])
  names(descendant_tips) <- internal_nodes
  descendant_tips
}

getDescTips(tree)



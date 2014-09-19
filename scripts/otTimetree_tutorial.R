## # Complete whale analysis workflow:

## ## Part 1. Load required libraries and functions

require(devtools)
#install_github("rotl", user="fmichonneau")
require(rotl)
require(rjson)
require(geiger)

## ## Load required input files
## The file CetaceanCompData.csv contains a list of taxa for which we have compartive data.
table <- read.csv("../data/CetaceanCompData.csv")
head(table)

## ## Download source trees from the Open Tree of Life database.
## Generate a list of taxa for which we have comparive data.
taxalist <- table$genspec
taxalist

## Query OpenTree api for taxonomic name resolution
tax <- rotl::tnrs_match_names(taxalist)
head(tax)

## Get ottIDs for the comparative data.
ottids <- as.character(tax$ott_id) 

## Find all source trees with our focal taxa in them, filter out trees without branch lengths. 
studyTrees <- findsourcetrees(ottids)      
possibleTrees <- getStudyTrees(studyTrees)
possibleTrees$summary

## We are going to choose study 1428 and tree2855. Under normal circumstances, you would choose 1927, because
## it is a timetree with all of our taxa. However, to test our method, we are going to use a more poorly matching
## tree, to demonstrate how we could match under non-ideal conditions.
tree <- extractTree(2, possibleTrees)
plot(tree)

## We'll also grab the other tree to validate our tree later.
reference_tree <- extractTree(1, possibleTrees)

## ## Now we will search the OT study database for trees calibrated in time units that shares taxa with our focal tree.
stcal <- sourceTreeCalibrations(tree, taxalist)
stcal$summary

## Get calibrations from other time calibrated timetrees. 
target <- drop.tip(tree, "Danio_rerio")
ot_calibrations <- getSourceTreeCalibrations(stcal)
head(ot_calibrations)

## Get fossil calibration data for our study group.
# Get the data (tree, occurrences)
names="cetacea"
pbdb.data<-pbdb_occurrences(limit="all",base_name=names,vocab="pbdb",show=c("phylo", "time", "ident"))
# Get genus + age related info
g <- getGenera(tree,pbdb.data)
# Assign calibrations to each node
fossil_calibrations<-getCalibrations(g)
# Removed nested calibrations
fixed.constraints<-filter.constraints(fossil_calibrations,tree)
# Assignning daughter nodes to calibrated parent nodes
calibrated.nodes<-fill.daughters(fixed.constraints,tree) 

## Combine calibrations into a single table
ot_calibrations$source <- rownames(ot_calibrations)
calibrated.nodes <- cbind(calibrated.nodes[,c(1,2,2, 5:6)])
calibrated.nodes$source <- "PBDB"
colnames(calibrated.nodes) <- colnames(ot_calibrations)
calibrations <- rbind(ot_calibrations, calibrated.nodes)
head(calibrations)

## Plot results
plot(target, cex=0.5)
nn <- sapply(1:nrow(calibrations), function(x) getMRCA(target, as.character(calibrations[x,4:5])))
nodelabels(node = nn, pch=21, bg=ifelse(calibrations$source=="PBDB", 3, 2), cex=log(calibrations$MaxAge)/3)

## Use pathd8 to get results. 
tr <- replaceHashes(target, calibrations)
phypd8 <- PATHd8.phylo(tr$scion, tr$table[1:23,], base = ".tmp_PATHd8", rm = FALSE)

## Plot time calibrated phylogeny
plot.phylo(phypd8, cex=0.5, tip.color = (phypd8$tip.label %in% taxalist)+1)
nodelabels(node = nn, pch=21, bg=ifelse(calibrations$source=="PBDB", 3, 2), cex=log(calibrations$MaxAge)/3)
#Vertical lines indicating 50 my intervals
ageMax <- max(nodeHeights(phypd8))
abline(v=ageMax - seq(0,ageMax, 50), lty=2)

## Compare final tree to results of known timetree. 
compTimeTree <- drop.tip(phypd8, phypd8$tip.label[!(phypd8$tip.label %in% taxalist)])
refTimeTree <- drop.tip(reference_tree, reference_tree$tip.label[!(reference_tree$tip.label %in% compTimeTree$tip.label)])
compTimeTree <- reorder(compTimeTree,"cladewise")
refTimeTree <- reorder(refTimeTree,"cladewise")
plot(c(0,0), c(0, 0), type="n", xlim=c(0, 60), ylim=c(0, 6), xaxt="n", yaxt="n",xlab="", ylab="", bty="n")
plotTree(refTimeTree, add=TRUE, color="red")
plotTree(compTimeTree, add=TRUE, ftype="off",lty=2, offset=c(max(nH2)-max(nH1),1))






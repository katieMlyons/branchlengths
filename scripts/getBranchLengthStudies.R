setwd("~/repos/OThackathon")
studydir <- "~/repos/phylesystem/shards/phylesystem-1/study/"

hook <- 'ot:substitutionCount'
outname <- "studies_substitionCount.txt"
system(paste("grep -r ", hook," ", studydir,"* >> ", outname, sep=""))
substCountStudies <- readLines(outname)
substCountTrees <- unname(sapply(substCountStudies, function(x) strsplit(x, ":")[[1]][1]))
writeLines(substCountTrees, con=outname)
substCountTaxa <- sapply(substCountTrees, function(y){x <- readLines(y); x[grep('ot:focalCladeOTTTaxonName', x)]})
tmp <-sapply(1:length(substCountTaxa), function(x) strsplit(substCountTaxa[[x]], ":")[1])
substCountTaxa <- sapply(tmp, function(x) x[length(x)])
substCountTaxa <- sapply(substCountTaxa, function(x) gsub("\"", "", x))
substCountTaxa <- sapply(substCountTaxa, function(x) gsub(" ", "", x))
substCountTaxa <- sapply(substCountTaxa, function(x) gsub(",", "", x))
substCountTaxaList <- unlist(substCountTaxa)


hook <- 'ot:changesCount'
outname <- "studies_changesCount.txt"
system(paste("grep -r ", hook," ", studydir,"* >> ", outname, sep=""))
changesCountStudies <- readLines(outname)
changesCountTrees <- unname(sapply(changesCountStudies, function(x) strsplit(x, ":")[[1]][1]))
writeLines(changesCountTrees, con=outname)
changesCountTaxa <- sapply(changesCountTrees, function(y){x <- readLines(y); x[grep('ot:focalCladeOTTTaxonName', x)]})
tmp <-sapply(1:length(changesCountTaxa), function(x) strsplit(changesCountTaxa[[x]], ":")[1])
changesCountTaxa <- sapply(tmp, function(x) x[length(x)])
changesCountTaxa <- sapply(changesCountTaxa, function(x) gsub("\"", "", x))
changesCountTaxa <- sapply(changesCountTaxa, function(x) gsub(" ", "", x))
changesCountTaxa <- sapply(changesCountTaxa, function(x) gsub(",", "", x))
changesCountTaxaList <- unlist(changesCountTaxa)

hook <- 'ot:time'
outname <- "studies_timeCount.txt"
system(paste("grep -r ", hook," ", studydir,"* >> ", outname, sep=""))
timeCountStudies <- readLines(outname)
timeCountTrees <- unname(sapply(timeCountStudies, function(x) strsplit(x, ":")[[1]][1]))
writeLines(timeCountTrees, con=outname)
timeCountTaxa <- sapply(timeCountTrees, function(y){x <- readLines(y); x[grep('ot:focalCladeOTTTaxonName', x)]})
tmp <-sapply(1:length(timeCountTaxa), function(x) strsplit(timeCountTaxa[[x]], ":")[1])
timeCountTaxa <- sapply(tmp, function(x) x[length(x)])
timeCountTaxa <- sapply(timeCountTaxa, function(x) gsub("\"", "", x))
timeCountTaxa <- sapply(timeCountTaxa, function(x) gsub(" ", "", x))
timeCountTaxa <- sapply(timeCountTaxa, function(x) gsub(",", "", x))
timeCountTaxaList <- unlist(timeCountTaxa)

hook <- '@length'
outname <- "studies_otherCount.txt"
system(paste("grep -r '", hook,"' ", studydir,"* >> ", outname, sep=""))
otherCountStudies <- readLines(outname)
otherCountTrees <- unname(sapply(otherCountStudies, function(x) strsplit(x, ":")[[1]][1]))
otherCountTrees <- unique(otherCountTrees)
allTrees <- unique(otherCountTrees)
otherCountTrees <- otherCountTrees[!(otherCountTrees %in% c(timeCountTrees, changesCountTrees, substCountTrees))]
writeLines(otherCountTrees, con=outname)
otherCountTaxa <- sapply(otherCountTrees, function(y){x <- readLines(y); x[grep('ot:focalCladeOTTTaxonName', x)]})
tmp <-sapply(1:length(otherCountTaxa), function(x) strsplit(otherCountTaxa[[x]], ":")[1])
otherCountTaxa <- sapply(tmp, function(x) x[length(x)])
otherCountTaxa <- sapply(otherCountTaxa, function(x) gsub("\"", "", x))
otherCountTaxa <- sapply(otherCountTaxa, function(x) gsub(" ", "", x))
otherCountTaxa <- sapply(otherCountTaxa, function(x) gsub(",", "", x))
otherCountTaxaList <- unlist(otherCountTaxa)


length(allTrees)

tmp <- readLines(otherCountTrees[626])
tmp[grep("ottTaxonName", tmp)]
otherCountTaxaList




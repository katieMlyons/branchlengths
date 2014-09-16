## Get a list
pantheria <- read.csv("./data/PantheriaReduced.csv")
Cetacea = subset(pantheria, Order=="Cetacea")
Cetacea <- subset(Cetacea, !is.na(Mass.g))
taxalist <- unique(Cetacea$genspec)
#Random sample for which we have data
samp <- sample(length(taxalist), 30, replace=FALSE)

datafile <- Cetacea[samp,]
rownames(datafile) <- NULL

write.csv(datafile, "./data/CetaceanCompData.csv")

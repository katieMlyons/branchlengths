
getGenera <- function(tree, pbdb.data) {
# this function will only work for trees containing tips of the format 
# Genus_species  
  
  genera<-data.frame(node=numeric(), genus=character(), age=numeric())
  for (i in (length(tree$tip.label)+1):(length(tree$tip.label)*2-1) ) {
      row=which(tree$edge[,1]==i)
      descendents=tree$edge[,2][row]
      dec1 <- getDescendants(tree, descendents[1])
      dec1 <- dec1[which(dec1 <= length(tree$tip.label))]
      dec2 <- getDescendants(tree, descendents[2])
      dec2 <- dec2[which(dec2 <= length(tree$tip.label))]
      
      s=strsplit(tree$tip.label[c(dec1, dec2)], "_")
      for(j in 1:length(s)) {
        genera<-rbind(genera,data.frame(node=i,genus=s[[j]][1],age=NA))
      }
    }
  
  for(k in as.character(unique(genera$genus))) {
    age <- pbdb.genus.min(pbdb.data, k)
    genera[which(k == genera$genus), "age"] <- age
  }
  
  return(unique(genera))
} # end of function
  
g <- getGenera(tree)

#######

# for each node, find out if it has more than one genus 
# listed more than one in table g
# if it does, take the oldest number from the age column
# 
# 

# assign that date, to the ancestral node
calibrations <- data.frame(node=unique(g$node), age=NA)
for(n in 1:nrow(calibrations)) {
  ages <- g[which(g$node == calibrations[n, "node"]), ]
  
  if(nrow(ages) > 0) {
    print(ages)
    max.age <- max(na.omit(ages$age))
    calibrations[n, "age"] <- max.age
  }

}



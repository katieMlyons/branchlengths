## returns unique tree ids and corresponding study ids for a set of ottIDs
findsourcetrees <- function(otts){
##start with a list of ottids


  findstudies <- function(otts){
    res <- rotl::studies_find_trees(property="ot:ottId", value=otts)
    res <- fromJSON(content(res, "text"))$matched_studies
  
  }

  gettrees <- function(x, studies){
    studyid <- studies[[x]]['ot:studyId']
    matched_trees <- trees[[x]]$matched_trees
    studs <-  lapply(1:length(matched_trees), function(x) cbind(studyid[x,], matched_trees[[x]]$nexson_id))
    df <- do.call(rbind, studs)
    colnames(df) <- c("studyid", "treeid")
    df
  }


studies <- lapply(otts, function(x) findstudies(x))
table <- do.call(rbind, lapply(1:length(studies), gettrees, studies=studies))
unique(table)

}
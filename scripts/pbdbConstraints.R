pbdb.fossil.min<-function(names) 

{
	
pbdb.data<-pbdb_occurrences(limit="all",base_name=names,vocab="pbdb",show=c("phylo", "time", "ident"))

if(!length(pbdb.data)==0) { 

early.ages=pbdb.data$early_age
max.age=max(early.ages)

}

else {
	max.age=0;
}

if(max.age!=0) {
	return(max.age) # this represents the minimum maximum age!
}
else {
	print("I'm sorry, no calibration info available!")

}
#EOF
}


### notes on installing & using the PBDB API in R can be found here:
### https://github.com/ropensci/paleobioDB

### Installing paleobioDB

### From CRAN:
### install.packages("paleobioDB")
### library(paleobioDB)

### From GitHub:
### install.packages("devtools")
### library(devtools)
### install_github("ropensci/paleobioDB")
### library(paleobioDB)

### paleobioDB functions are (non-exhaustively) described here: http://cran.r-project.org/web/packages/paleobioDB/paleobioDB.pdf
### Furhter details can be found here: http://paleobiodb.org/data1.1/occs/list

### Using the PBDB API to download the occurrence data: an example using Cetacea

### A note from PBDB:
### Beware of synonyms and errors, they could twist your estimations about species richness, evolutionary and extinction rates, etc. paleobioDB users should be critical about the raw data downloaded from the database and filter the data before analyzing it.

### To download occurrence data for a given group:
	
### cetacea<-pbdb_occurrences(limit="all",base_name="cetacea",vocab="pbdb",show=c("phylo", "time", "ident"))

### Options for specifying taxa
### taxon_name: Return only records associated with the specified taxonomic name(s). You may specify multiple names, separated by commas.
### base_name: Return only records associated with the specified taxonomic name(s), *or any of their children*. 
### You may specify multiple names, separated by commas. -> THIS OPTION DOESN'T SEEM TO WORK -> ACTUALLY, I THINK IT DOES
### Note PBDB does not generally handle common names e.g. use Canis, not dog

### Options for headers
### vocab="pbdb" : The original Paleobiology Database field names. This vocabulary is the default for text format responses (.tsv, .csv, .txt).
### vocab="com" : 3-character abbreviated field names. This is the default for JSON responses.

### Options for time periods
### interval: Return only records whose temporal locality falls within the named geologic time interval.

### show: specify additional fields associated with the occurrence/collection
### coords: The latitude and longitude of an occurrence
### phylo: Additional information about the taxonomic classification of the occurence
### ident: The actual taxonomic name by which this occurrence was identified
### strat: Basic information about the stratigraphic context of the occurrence
### time: Additional information about the temporal locality of the occurrence
	### early_age: The early bound of the geologic time range associated with this occurrence (in Ma)
	### late_age: The late bound of the geologic time range associated with this occurrence (in Ma)

### We need info. about the (1) the taxonomy, & (2) the stratigraphy.
### show "ident" required to include genes_name
### show "time" required to include early_age

### To view the headers

### names(cetacea)

### To obtain a list of genera species, or minimum age estimates

### cetacea$genus_name
### cetacea$species_name
### cetacea$late_age

### [given two sets of exchangability taxa...]

### 1. define the family (if available) or a list of exchangeable taxa
### names="Balaenidae"

### names=c("Balaenidae","Eubalaena","Idiocetus","Morenocetus","Palaeocetus","Balaena","Balaenotus","Balaenula","Balaenella","Probalaena") # this should give you the same result

### pbdb.data<-pbdb_occurrences(limit="all",base_name=names,vocab="pbdb",show=c("phylo", "time", "ident"))

### 2a. find the oldest minimum age
### late_ages=pbdb.data$late_age
### max.age=max(late_age)

# ----

### 2b. find the oldest minimum age for a sub group
### row=which(cetacea$genus_name=="Balaena")
### late_ages=pbdb.data$late_age[row]
### max.age=max(late_ages)


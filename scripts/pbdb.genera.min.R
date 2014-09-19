


pbdb.genus.min<-function(pbdb.data,genera,upper.bound=100)
# This function gets the youngest secure age of the oldest fossil occurrence for a given genera.
# A default is specified for the upper bound but you may with to change this depending on the group.
  
{

  if(!length(pbdb.data)==0) { # if the data frame is not empty
    
    # extract columns for genus only
    rows=which(pbdb.data$genus_name==genera)    
    late.ages=pbdb.data$late_age[rows] # colunm contains: late bound of the geologic time range associated with this occurrence (in Ma)
    late.ages=sort(late.ages) # sort a -> b
    
    ### the following is a quick solution
    ### for eliminating dates that are implausibly old
    count=0
    
    while(count==0) {
      
      if (length(late.ages) > 0) { # if there are available dates
        
        max.age=max(late.ages) # define the greatest
        
        if(max.age >=upper.bound) { # if this is greater than the specified upper limit
          late.ages<-late.ages[-length(late.ages)] # remove the last element of this list
        }
        else {
          
          count=1	
        }
      }
      # this would mean that there were no suitable dates available in the list
      else {
        max.age=0
        count=1
      }	
    }
    ###
  }
  
  else {
    max.age=0;
  }
  
  if(max.age!=0) {
    return(max.age) # this represents the minimum maximum age!
  }
  else {
    return(NA)
#     print("I'm sorry, no calibration info available!")
    
  }
  #EOF
}

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")

library(data.table)
library(dplyr) 

stations <- fread(paste0(getwd(),"/", "available_stations_3.csv"))

################################################################################
# example with one station for testing
exa <- stations[2,]$ID
exa <- fread(paste0(getwd(),"/", exa,".csv"))
exa <- exa %>% filter(element %in% c("TMAX","TMIN","PRCP"))

means <- exa %>% 
         select(-c(1,2)) %>% 
         group_by(element,month) %>% 
         summarise(across(everything(), .f = list(mean = mean), na.rm = TRUE))       

means <- data.frame(cbind(means[,c(1,2)],lapply(means[,c(-1,-2)], function(x) round(x, digits = 0))))

# unfortunately there are some means that are NA, so we need to replace that NA
for (i in 1:(12*3)){
  pos <- which(!is.finite(unlist(means[i,c(-1,-2)])))
  for (j in pos){
    if(j != 1 ){
      means[i,2+j] <- means[i,j + 1]
    }
    else if(i %in% c(1,13,25)){
      means[i,2+j] <- means[i - 1,31 + 2]
    }else{
      means[i,2+j] <- means[i + 11,31 + 2]
    }   
  }
}

check <- unlist(lapply(means[,c(-1,-2)], function(x) all(is.finite(x))))
if(!all(check)){
  print("there still are NAs")
}

pos_na <- sapply(exa[,c(-1,-2,-3,-4)], function(x) which(!is.finite(x) == "TRUE"))

for(i in 1:31){ 
  a <- left_join(exa[unlist(pos_na[i]),], means, by=c('element','month')) %>% 
                                    select(paste("Val",i,"_mean",sep="")) %>% 
                                               arrange(.by_group = FALSE) %>%
                                                   as.vector()
  b <- unlist(pos_na[i])
  if(length(b) > 0){
    for(j in 1:length(b)){
      exa[b[j]][[(i+4)]] <- as.numeric(a[j])
    }
  }
}

pos_na_2 <- sapply(exa[,c(-1,-2,-3,-4)], function(x) which(!is.finite(x) == "TRUE"))




# this function computes the mean for each variable(TMAX,TMIN,PRCP), each month and each day.
# Then it adds the mean to the right missing values considering the previous three positions
# in the table.
# It has problems if all the means in a month and all the previous are NAs  
fillMissingValues <- function(station_ID){  # input must be a string
  exa <- fread(paste0(getwd(),"/", station_ID,".csv"))
  exa <- exa %>% filter(element %in% c("TMAX","TMIN","PRCP"))
  
  means <- exa %>% 
    select(-c(1,2)) %>% 
    group_by(element,month) %>% 
    summarise(across(everything(), .f = list(mean = mean), na.rm = TRUE))       
  
  means <- data.frame(cbind(means[,c(1,2)],lapply(means[,c(-1,-2)], function(x) round(x, digits = 0))))
  for (i in 1:(12*3)){
    pos <- which(!is.finite(unlist(means[i,c(-1,-2)])))
    for (j in pos){
      if(j != 1 ){
        means[i,2+j] <- means[i,j + 1]
      }
      else if(i %in% c(1,13,25)){
        means[i,2+j] <- means[i - 1,31 + 2]
      }else{
        means[i,2+j] <- means[i + 11,31 + 2]
      }   
    }
  }
  
  check <- unlist(lapply(means[,c(-1,-2)], function(x) all(is.finite(x))))
  if(!all(check)){
     print("there still are NAs among the means. So, watch out")
  }
  
  pos_na <- sapply(exa[,c(-1,-2,-3,-4)], function(x) which(!is.finite(x) == "TRUE"))
  
  for(i in 1:31){ 
    a <- left_join(exa[unlist(pos_na[i]),], means, by=c('element','month')) %>% 
      select(paste("Val",i,"_mean",sep="")) %>% 
      arrange(.by_group = FALSE) %>%
      as.vector()
    b <- unlist(pos_na[i])
    if(length(b) > 0){
      for(j in 1:length(b)){
        exa[b[j]][[(i+4)]] <- as.numeric(a[j])
      }
    }
  }
  
  pos_na_2 <- sapply(exa[,c(-1,-2,-3,-4)], function(x) which(!is.finite(x) == "TRUE"))
  
  summ = 0
  for (i in 1:31){ summ = summ + length(unlist(pos_na_2[i]))}
  if(summ != 0){
    print("there may be some problems. Number of missing values: " + summ)
  }
  fwrite(exa,  paste0(getwd(), "/1_",station_ID,".csv"))
  exa
}

#a <- fillMissingValues("USC00148235")
#mam2 <- fread(paste0(getwd(),"/1_", "USC00148235",".csv"))

for (i in stations$ID){
  print(i)
  fillMissingValues(i)
}



# useful file to obtain the stations for Texas given the list from Eno

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/")

library(data.table)
library(dplyr) 
library(tidyr)

stnscsv <- paste0(getwd(),"/../","stations.csv")

typedcols <- c( "A11", "F9", "F10", "F7", "X1","A2",
                "X1","A30", "X1", "A3", "X1", "A3", "X1", "A5" )
stns <- read.fortran(paste0(getwd(),"/","ghcnd-stations.txt"),
                     typedcols, 
                     comment.char="")

# Import the stations to analyze from Eno's results
final_tmax <- fread("final_tmax.csv")
final_tmin <- fread("final_tmin.csv")
final_prec <- fread("final_prcp.csv")

# join to add lat,long, alt and name for each station
final_tmax <- left_join(final_tmax, stns, by=c('ID'))[,c(2,5,6,7,9)]
final_tmin <- left_join(final_tmin, stns, by=c('ID'))[,c(2,5,6,7,9)]
final_prec <- left_join(final_prec, stns, by=c('ID'))[,c(2,5,6,7,9)]

#################################################################################
# i need to extract the file with these stations for texas 

list_stations <- unique(c(final_tmax$ID,final_tmin$ID,final_prec$ID))

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project")

library(data.table)
library(dplyr) 

numFiles <- length(list_stations)
dirname <- paste0(getwd(),"/","ghcnd_all/ghcnd_all/")

for (i in 1:numFiles) {
  
  infile <- paste0(dirname, list_stations[i], ".dly")
  outfile <- paste0(getwd(),"/experiment/texas/", list_stations[i], ".csv")
  
  cols <- c( "A11", "I4", "I2", "A4",
             rep( c( "I5", "A1", "A1", "A1"), 31) )
  df <- read.fortran(infile, cols, na.strings="-9999") # -9999 indicates missing data
  
  # next, fill in the column names
  tmp <- c("Val","xxM","xxQ","xxS") # xx so we can ditch them later
  vhdrs <- paste(   rep(tmp,31),   rep(1:31,each=4), sep="")
  hdrs <- c("ID", "year", "month", "element", vhdrs)
  names(df) <- hdrs
  df <- df[df$year >= 1965 & df$year <= 2015,]
  df_out <- dplyr::select(df, -matches("xx*")) # get rid of M, Q, S 
  fwrite(df_out, outfile)
}

#################################################################################
# now i need to fill the missing values using missingValues.R
# notice that the stations without missing values are saved with a new name 1_stationId

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/texas")
#list_stations <- paste0("1_",list_stations)

list <- vector()

for (i in list_stations){
  print(i)
  list[i] <- fillMissingValues(i)
}
#################################################################################




########
#######  still woriking on the following lines
########






# Create a daily dataframe of 1 year of a feature (TMAX,TMIN,PRCP) as mean over 1965-2015

df_final1 = df %>% select(month,day,daycount)

numFiles <- dim(stns)[1]

for (i in 1:numFiles) { # For each station
  fname <- paste0(getwd(),"1_",stns$ID[i],"_50.csv")
  tempf <- fread(fname, stringsAsFactors=FALSE)
  # tempf = tempf[tempf$year>=1965 & tempf$year<=2015,]     # Select range of interest (default 1965-2015) 
  df = tempf %>% group_by(month,day) %>% summarise(feat = mean(PRCP, na.rm=TRUE))
  df = df[complete.cases(df),]                      # Removes non existing dates (eg. 31 Feb), now there are 366 days
  colnames(df) = c("month","day",first(tempf$ID))   # Creates a day counter from 1 to 366
  df$daycount = seq.int(nrow(df))
  
  df_final1 = merge(df_final1,df,by = c("month","day","daycount"))
  
}

df_final1 = df_final1 %>% arrange(daycount) # Order datatset with respect to daycount (cronological order)

outfile1 <- paste0(getwd(),"PRCP_tot_mean.csv")
fwrite(df_final1, outfile1)

##################################################################################################################
##################################################################################################################

# Create monthly dataset of a feature as a mean over 1965-2015
numFiles <- dim(stns)[1]

df_final2 = data.frame(matrix(nrow = 12,ncol= 0))

for (i in 1:numFiles) { # For each station
  fname <- paste0("C:/Users/enogj/Desktop/proj stat/dataset/Ordered dataset/",stns$ID[i],"_50.csv")
  tempf <- fread(fname, stringsAsFactors=FALSE)
  # tempf = tempf[tempf$year>=1965 & tempf$year<=2015,]     # Select range of interest (default 1965-2015)
  df = tempf %>% group_by(month) %>% summarise(feat = mean(PRCP, na.rm=TRUE))
  df = df %>% select(feat)
  colnames(df) = c(first(tempf$ID))
  df_final2 = cbind(df_final2,df)
  
}

outfile1 <- paste0("C:/Users/enogj/Desktop/proj stat/dataset/PRCP/","PRCP_month_mean.csv")
fwrite(df_final2, outfile1)

##################################################################################################################
##################################################################################################################

# Creates a annual dataset for kriging as mean over 1965-2015

numFiles <- dim(stns)[1]

df_final3 = data.frame()


for (i in 1:numFiles) { # For each station
  fname <- paste0("C:/Users/enogj/Desktop/proj stat/dataset/Ordered dataset/",stns$ID[i],"_50.csv")
  tempf <- fread(fname, stringsAsFactors=FALSE)
  # tempf = tempf[tempf$year>=1965 & tempf$year<=2015,]     # Select range of interest (default 1965-2015)
  df <- tempf %>% group_by(ID) %>% summarise(ID = first(ID),value = mean(TMAX, na.rm=TRUE))
  df$LAT=stns$LAT[i]
  df$LON=stns$LON[i]
  df$ELEV=stns$ELEV[i]
  df_final3 = rbind(df_final3,df)
  
}

outfile1 <- paste0("C:/Users/enogj/Desktop/proj stat/dataset/TMAX/","TMAX_tot_mean_krg.csv")
fwrite(df_final3, outfile1)

##################################################################################################################
##################################################################################################################

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/")

# useful file to obtain the stations for Texas given the list from Eno
library(data.table)
library(dplyr) 
library(tidyr)

stnscsv <- paste0(getwd(),"stations.csv")

typedcols <- c( "A11", "F9", "F10", "F7", "X1","A2",
                "X1","A30", "X1", "A3", "X1", "A3", "X1", "A5" )
stns <- read.fortran(paste0(getwd(),"/","ghcnd-stations.txt"),
                     typedcols, 
                     comment.char="")
hdrs <- c("ID", "LAT", "LON", "ELEV", "ST", "NAME","GSN", "HCN", "WMOID")
names(stns) <- hdrs
setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/")

# Import the stations to analyze from Eno's results
final_tmax <- fread("final_tmax.csv")
final_tmin <- fread("final_tmin.csv")
final_prcp <- fread("final_prcp.csv")

# join to add lat,long, alt and name for each station
final_tmax <- left_join(final_tmax, stns, by=c('ID'))[,c(2,5,6,7,9)]
final_tmin <- left_join(final_tmin, stns, by=c('ID'))[,c(2,5,6,7,9)]
final_prcp <- left_join(final_prcp, stns, by=c('ID'))[,c(2,5,6,7,9)]

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/texas/")
# recap tables
fwrite(final_tmax,"final_tmax_f.csv")
fwrite(final_tmin,"final_tmin_f.csv")
fwrite(final_prcp,"final_prcp_f.csv")

#################################################################################
# i need to extract the file with these stations for texas 

# Import the stations to analyze from Eno's results
final_tmax <- fread("final_tmax_f.csv")
final_tmin <- fread("final_tmin_f.csv")
final_prcp <- fread("final_prcp_f.csv")

list_stations <- unique(c(final_tmax$ID,final_tmin$ID,final_prcp$ID))

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


setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/texas/")

# function to create the table as we would like to have
getTable <- function(est_states,string){ #string must be one among TMAX, TMIN and PRCP
  
  numFiles <- dim(est_states)[1]
  for (i in 1:numFiles) {

    df <- data.frame( year = rep(c(1965:2015), each = 12), month = rep(c(1:12),51))
    df$ID <- est_states$ID[i]
    
    fname <- paste0(getwd(),"/1_",est_states$ID[i],".csv")
    stnf <- fread(fname)  ## read the file for the current station
   
    ## select only the TMAX elements from the incoming file, and merge with the data frame df.
    ## Use gather (from tidyr) to convert the daily data from row to columns.
    outfr <- stnf[stnf$element == string, -c(4)]
    outfr <- merge(df, outfr, all.x=TRUE)
    outfr1 <- outfr %>% gather(day, TMAX, Val1:Val31)  ## create one line for each day
    outfr1$day <- as.integer(substring(outfr1$day, 4, last = 1000000L)) 
    ## label the day by the integer after "Val". last = 1000000L means go to the end of the string.
    outfr1 <- outfr1 %>% arrange(year,month,day)
    
    ## save the file
    outfile <- paste0( getwd(),"/",string,"/2_",est_states$ID[i],".csv")
    fwrite(outfr1, outfile)
  }
}

getTable(final_tmax,"TMAX")
getTable(final_tmin,"TMIN")
getTable(final_prcp,"PRCP")


################################################################################
################################################################################









# Create a daily dataframe of 1 year of a feature (TMAX,TMIN,PRCP) as mean over 1965-2015


numFiles <- dim(stns)[1]

for (i in 1:numFiles) { # For each station
  fname <- paste0(getwd(),"/1_",stns$ID[i],".csv")
  tempf <- fread(fname, stringsAsFactors=FALSE)
  # tempf = tempf[tempf$year>=1965 & tempf$year<=2015,]     # Select range of interest (default 1965-2015) 
  df = tempf %>% group_by(month,day) %>% summarise(feat = mean(PRCP, na.rm=TRUE))
  df = df[complete.cases(df),]                      # Removes non existing dates (eg. 31 Feb), now there are 366 days
  colnames(df) = c("month","day",first(tempf$ID))   # Creates a day counter from 1 to 366
  df$daycount = seq.int(nrow(df))
  
  df_final1 = merge(df_final1,df,by = c("month","day","daycount"))
  
}

df_final1 = df %>% select(month,day,daycount)

df_final1 = df_final1 %>% arrange(daycount) # Order datatset with respect to daycount (cronological order)

outfile1 <- paste0(getwd(),"PRCP_tot_mean.csv")
fwrite(df_final1, outfile1)

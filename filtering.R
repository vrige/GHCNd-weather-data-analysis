setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project")

########################################################################
#
# Extraction dataset stations.csv and inventory.csv
# Plus selecting 2,000 stations in the contiguous US that have data
# spanning the 80 years from 1936-2015.
#
########################################################################


stnscsv <- paste0(getwd(),"/","stations.csv")

typedcols <- c( "A11", "F9", "F10", "F7", "X1","A2",
                "X1","A30", "X1", "A3", "X1", "A3", "X1", "A5" )
stns <- read.fortran(paste0(getwd(),"/","ghcnd-stations.txt"),
                     typedcols, 
                     comment.char="")
hdrs <- c("ID", "LAT", "LON", "ELEV", "ST", "NAME","GSN", "HCN", "WMOID")
names(stns) <- hdrs
write.csv(stns,stnscsv)

head(stns)

########################################################################
########################################################################

inventorycsv <- paste0(getwd(),"/","inventory.csv")

invcols <- c( "A11", "X1", "F8", "X1", "F9", "X1","A4",
              "X1","I4", "X1", "I4" )
inv <- read.fortran(paste0(getwd(),"/","ghcnd-inventory.txt"),
                    invcols,
                    comment.char="")
invhdrs <- c("ID", "LAT", "LON", "ELEM" , "FIRST", "LAST")
names(inv) <- invhdrs
write.csv(inv,inventorycsv)

head(inv)

########################################################################
########################################################################
library(dplyr)


#select all stations in the US 
us1 <- inv[grep("US+",inv$ID),]    #maybe it's better "^US"

us1 <- us1 %>% mutate(period = LAST - FIRST + 1)

us1_stat <- us1 %>%
  summarise(first = mean(us1$FIRST),LAST = mean(us1$LAST),
            mean_working_time = mean(period), max = max(period), n = n())

us1_stat2 <- us1 %>%
  filter((period) > us1_stat$mean_working_time)

#plot(us1$LON,us1$LAT)

# stations that monitor max temperature readings (TMAX)
us2 <- us1[us1$ELEM == "TMAX",]
# between 1936 and 2015
us_years <- us2[(us2$FIRST <= 1965) & us2$LAST >= 2015, c(1,5,6)]
colnames(us_years) <- c("ID", "TMAXf", "TMAXl")

# and have also been monitoring TMIN between 1936 and 2015
us3 <- us1[(us1$ELEM == "TMIN") & (us1$FIRST <= 1965) & (us1$LAST >= 2015), c(1,5,6)]
colnames(us3) <- c("ID", "TMINf", "TMINl")
us_years <- merge(us_years, us3, all=FALSE)

# in addition to monitoring PRCP between 1936 and 2015
us4 <- us1[(us1$ELEM == "PRCP") & (us1$FIRST <= 1965) & (us1$LAST >= 2015), c(1,5,6)]
colnames(us4) <- c("ID", "PRCPf", "PRCPl")
us_years <- merge(us_years, us4, all=FALSE)

us5 <- us1[(us1$ELEM == "SNOW") & (us1$FIRST <= 1965) & (us1$LAST >= 2015), c(1,5,6)]
colnames(us5) <- c("ID", "SNOWf", "SNOWl")
us_years <- merge(us_years, us5, all=FALSE)

us6 <- us1[(us1$ELEM == "SNWD") & (us1$FIRST <= 1965) & (us1$LAST >= 2015), c(1,5,6)]
colnames(us6) <- c("ID", "SNWDf", "SNWDl")
us_years <- merge(us_years, us6, all=FALSE)

us_stns <- stns[grep("US+",stns$ID),]
unique(us_stns$ST)

# [1] "SD" "CO" "NE" "AK" "AL" "AR" "AZ" "CA" "TN" "CT" "DC" "DE" "FL" "GA" "HI" "IA" "ID" "IL" "IN" "KS"
#[21] "KY" "LA" "MA" "MD" "ME" "MI" "MN" "MO" "MS" "MT" "NC" "ND" "NH" "NJ" "NM" "NV" "NY" "OH" "OK" "OR"
#[41] "PA" "RI" "SC" "TX" "UT" "VA" "VT" "WA" "WI" "WV" "WY" "PI" "UM"

# PI stands for Pacific Islands, UM for US Minor Outlying Islands

conus_stns <-  us_stns[(us_stns != "AK") & 
                         (us_stns != "HI") & 
                         (us_stns != "PI") & 
                         (us_stns != "UM"), ]
us80 <- merge(conus_stns, us_years, all.y=TRUE)


########################################################################
########################################################################

library(data.table)
library(dplyr) 
numFiles <- length(us80$ID)
dirname <- paste0(getwd(),"/","ghcnd_all/ghcnd_all/")

for (i in 1:numFiles) {
  
  infile <- paste0(dirname, us80$ID[i], ".dly")
  outfile <- paste0(getwd(),"/experiment/", us80$ID[i], ".csv")
  
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

########################################################################
########################################################################

us80$natmax <- 0
us80$natmin <- 0
us80$naprcp <- 0
us80$nasnow <- 0
us80$nasnwd <- 0
numFiles <- length(us80$ID)

for (i in 1:numFiles) {
  infile <-  paste0(getwd(),"/experiment/", us80$ID[i], ".csv")
  df <- read.csv(infile)
  us80$natmax[i] <- sum(sapply(df[df$element == "TMAX",c(5:35)], function (x) length(x[is.na(x)])))
  us80$natmin[i] <- sum(sapply(df[df$element == "TMIN",c(5:35)], function (x) length(x[is.na(x)])))
  us80$naprcp[i] <- sum(sapply(df[df$element == "PRCP",c(5:35)], function (x) length(x[is.na(x)])))
  us80$nasnow[i] <- sum(sapply(df[df$element == "SNOW",c(5:35)], function (x) length(x[is.na(x)])))
  us80$nasnwd[i] <- sum(sapply(df[df$element == "SNWD",c(5:35)], function (x) length(x[is.na(x)])))
}

# 2% of data = 50 x 12 x 31 x 0.05 = 372
thr = 50 * 12 * 31 * 0.02
stn80 <- us80[us80$natmax <= thr & us80$natmin <= thr & us80$naprcp <= thr & us80$nasnow <= thr & us80$nasnwd <= thr, ]
fwrite(stn80,  paste0(getwd(),"/experiment/", "stn80.csv"))

dim(stn80)
setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")




#######################################################################################
#
setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")

library(data.table)
library(dplyr) 

dirname <- paste0(getwd(),"/stn80.csv")
stn80 <- fread(dirname)
numFiles <- dim(stn80)[1]
df <- data.frame()

for (i in 1:numFiles) {
  
  print(i)
  infile <- paste0(getwd(),"/", stn80$ID[i], ".csv")
  temp <- fread(infile)
  df <- rbind(df,temp)
  
}
df[is.na(df)] <- -9999
fwrite(df,  paste0(getwd(), "/stn80_all.csv"))


################################################################################

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")

library(data.table)
library(dplyr) 

dirname <- paste0(getwd(),"/stn80.csv")
stn80 <- fread(dirname)
numFiles <- dim(stn80)[1]
df_tmax <- data.frame()
df_tmin <- data.frame()
df_prcp <- data.frame()
df_snow <- data.frame()
df_snwd <- data.frame()

for (i in 1:numFiles) {

  print(i)
  infile <- paste0(getwd(),"/", stn80$ID[i], ".csv")
  temp <- fread(infile)
  #temp[is.na(temp)] <- -9999
  tmax <- temp %>% filter(element == "TMAX")
  tmin <- temp %>% filter(element == "TMIN")
  prcp <- temp %>% filter(element == "PRCP")
  snow <- temp %>% filter(element == "SNOW")
  snwd <- temp %>% filter(element == "SNWD")
  df_tmax <- rbind(df_tmax,tmax)
  df_tmin <- rbind(df_tmin,tmin)
  df_prcp <- rbind(df_prcp,prcp)
  df_snow <- rbind(df_snow,snow)
  df_snwd <- rbind(df_snwd,snwd)
  
}

fwrite(df_tmax,  paste0(getwd(), "/stn80_tmax_na.csv"))
fwrite(df_tmin,  paste0(getwd(), "/stn80_tmin_na.csv"))
fwrite(df_prcp,  paste0(getwd(), "/stn80_prcp_na.csv"))
fwrite(df_snow,  paste0(getwd(), "/stn80_snow_na.csv"))
fwrite(df_snwd,  paste0(getwd(), "/stn80_snwd_na.csv"))


################################################################################
setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")

library(data.table)
library(dplyr) 
library(tidyr)

dirname <- paste0(getwd(),"/stn80.csv")
stn80 <- fread(dirname)
numFiles <- dim(stn80)[1]

#TMAX

infile <- paste0(getwd(),"/stn80_tmax_na.csv")
df <- fread(infile)
  
a <- df %>% group_by(ID,year,month) %>% summarise(across(c(2:32), ~ sum(is.na(.x))))
a$na_month <- rowSums(a[4:34])
a <- a[,c(1,2,3,35)]
b <- a %>% spread(month,na_month)
b[is.na(b)] <- 31
b$na_year <- rowSums(b[3:14])
b <- b[,c(1,2,15)]
c <- b %>% group_by(ID) %>% 
         spread(year, na_year, fill = 372) %>% 
         gather(2:52, key = "year", value = "na_year") %>%
         arrange(ID) 
d <- c %>% group_by(ID) %>% 
  spread(year, na_year) %>% 
  ungroup
d$na_id <- rowSums(d[2:52])
d <- d[,c(1,53)]
d$est_missing_years <- d[,2]/372
su <- data_frame(missing_year_threshold=integer(),n = integer(),fraction = double())
for (year in 1:50) { 
  su <- rbind(su,d %>% filter(est_missing_years > year) %>%  summarise(missing_year_threshold = year, n = n(), fraction = n()/274))
}



# this function compute the number of missing values for each station
# the difference with the tutorial is that we count also the period not mentioned
# in the data. (e.g. if 1 year is missing, we add 372 missing values)
# Modifying the output of the function is possible to extract the exact location of
# missing values by years and months
missing_values <- function(infile) {
  df <- fread(infile)
  a <- df %>% group_by(ID,year,month) %>% 
    summarise(across(c(2:32), ~ sum(is.na(.x))))
  a$na_month <- rowSums(a[4:34])
  a <- a[,c(1,2,3,35)]
  b <- a %>% spread(month,na_month)
  b[is.na(b)] <- 31
  b$na_year <- rowSums(b[3:14])
  b <- b[,c(1,2,15)]
  c <- b %>% group_by(ID) %>% 
    spread(year, na_year, fill = 372) %>% 
    gather(2:52, key = "year", value = "na_year") %>%
    arrange(ID) 
  d <- c %>% group_by(ID) %>% 
    spread(year, na_year) %>% 
    ungroup
  d$na_id <- rowSums(d[2:52])
  d <- d[,c(1,53)]
  d$est_missing_years <- d[,2]/372
  return(d)
}


infile <- paste0(getwd(),"/stn80_tmax_na.csv")
tmax <- missing_values(infile)
tmax_1y <- tmax %>% filter(est_missing_years < 1)
colnames(tmax_1y) <- c("ID", "tmax_na_id", "est_missing_years")


infile <- paste0(getwd(),"/stn80_tmin_na.csv")
tmin <- missing_values(infile)
tmin_1y <- tmin %>% filter(est_missing_years < 1)
colnames(tmin_1y) <- c("ID", "tmin_na_id", "est_missing_years")


infile <- paste0(getwd(),"/stn80_prcp_na.csv")
prcp <- missing_values(infile)
prcp_1y <- prcp %>% filter(est_missing_years < 1)
colnames(prcp_1y) <- c("ID", "prcp_na_id", "est_missing_years")


infile <- paste0(getwd(),"/stn80_snow_na.csv")
snow <- missing_values(infile)
snow_1y <- snow %>% filter(est_missing_years < 1)
colnames(snow_1y) <- c("ID", "snow_na_id", "est_missing_years")


infile <- paste0(getwd(),"/stn80_snwd_na.csv")
snwd <- missing_values(infile)
snwd_1y <- snwd %>% filter(est_missing_years < 1)
colnames(snwd_1y) <- c("ID", "snwd_na_id", "est_missing_years")


available_stations <- merge(tmax_1y[,1:2], tmin_1y[,1:2], all=FALSE)
available_stations <- merge(available_stations, prcp_1y[,1:2], all=FALSE)
available_stations2 <- merge(available_stations, snow_1y[,1:2], all=FALSE)
available_stations2 <- merge(available_stations2, snwd_1y[,1:2], all=FALSE)

fwrite(available_stations,  paste0(getwd(),"/", "available_stations_3.csv"))
fwrite(available_stations2,  paste0(getwd(),"/", "available_stations_5.csv"))



################################################################################
################################################################################
# some corrections and graphics of the considered area

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment")

library(data.table)
library(dplyr) 
library(tidyr)

dirname <- paste0(getwd(),"/stn80.csv")
stn80 <- fread(dirname)
stations <- fread(paste0(getwd(),"/", "available_stations_3.csv"))

stations_ext <- inner_join(stations,stn80,by="ID")[,c(1,5,6,7,8)]

stations <- stations_ext[((stations_ext$ST != "AK") & 
                            (stations_ext$ST != "HI") & 
                            (stations_ext$ST != "PI") & 
                            (stations_ext$ST != "UM")), ]

#fwrite(stations,  paste0(getwd(),"/", "available_stations_3.csv"))

library(ggplot2)
library(maps)
library(mapproj)

head(map_data('state'),15)
states <- c('OH','IN','MI','IL','WI','KY','WV')
est_states = stations[ (stations$ST %in% states) & (stations$ID != 'USW00003859'),]

usa_tbl <- map_data("state", region = c('ohio','indiana','michigan:south','illinois','wisconsin')) %>% as_tibble()

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "gray30", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = est_states, aes(x=LON.y,y=LAT.y), color = "red")












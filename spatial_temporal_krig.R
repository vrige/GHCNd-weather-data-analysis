library(gstat)
library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos) 
library(dplyr)

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/texas/")

final_tmax <- fread("final_tmax_f.csv")
final_tmin <- fread("final_tmin_f.csv")
final_prcp <- fread("final_prcp_f.csv")

# we want to save the variable time in UNIX format
# it starts from 1st January 1970, so i will discard the data from 1965 to 1970

day <- 60 *60 * 24 # seconds in a day
days <- 372 * 46 # days in 45 years (notice that all the months have 31 days)
end <- day  * days
tempo <- seq(0, end, by = day) 
tempo2 <- as.POSIXlt(tempo, origin="1970-01-01")
tempo2 <- as.character.POSIXt(tempo2[1:dim])

stns <- final_tmax$ID
data <- data.frame()

for (i in stns){
  temp <- fread(paste0(getwd(),"/TMAX/2_",i,".csv"))
  temp <- temp %>% filter(year >= 1970)
  
  # watch out the following line
  temp <- temp %>% filter(year == 1970)
  
  dim <- length(temp$year)
  temp$time <- tempo2[1:dim]
  
  
  data <- rbind(data,temp)
}
data <- data %>% dplyr::select(ID,TMAX,time)
data %>% filter(time == 0)

data <- left_join(data,final_tmax,by=c("ID"))
data <- data %>% dplyr::select(ID,TMAX,time,LAT,LON,ELEV)

#coordinates(data) <- c('LON','LAT')#,'ELEV')






#Create a SpatialPointsDataFrame
coordinates(data)=~LON+LAT
projection(data)=CRS("+init=epsg:4326")

#Transform into Mercator Projection
data.UTM <- spTransform(data,CRS("+init=epsg:3395")) 

dataSP <- SpatialPoints(data.UTM@coords,CRS("+init=epsg:3395")) 
#dupl <- zerodist(dataSP) 

dataDF <- data.frame(tmax=data.UTM$TMAX) 

dataTM <- as.POSIXct(data.UTM$time,tz="CET")

timeDF <- STIDF(dataSP,dataTM,data=dataDF) 

x11()
stplot(timeDF) 

var <- variogramST(tmax~1,data=timeDF,tunit="days",assumeRegular=F)#,na.omit=T) 




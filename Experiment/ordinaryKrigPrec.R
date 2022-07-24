## Exploratory Analysis ##
##----------------------##
library(data.table)

setwd("C:/Users/39340/Desktop/poliMI/Applied statistics/project/experiment/eno_colors")

sea=data.frame(LAT=27.5, LON=-96) 
cord.dec = SpatialPoints(cbind(sea$LON, sea$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))

cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
cord.UTM = as.data.frame(cord.UTM)

sea$x = cord.UTM$coords.x1
sea$y = cord.UTM$coords.x2

coordinates(sea)=c('x','y')

prcp_2<- fread(paste0(getwd(),"/PRCP/PRCP_07-13_mean_krg.csv"))
prcp_1<- fread(paste0(getwd(),"/PRCP/PRCP_87-93_mean_krg.csv"))
prcp_0<- fread(paste0(getwd(),"/PRCP/PRCP_67-73_mean_krg.csv"))

prcp_0 <- prcp_0 %>% filter(prcp_0$ELEV >= 0)
prcp_1 <- prcp_1 %>% filter(prcp_1$ELEV >= 0)
prcp_2 <- prcp_2 %>% filter(prcp_2$ELEV >= 0)

cord.dec_0 = SpatialPoints(cbind(prcp_0$LON, prcp_0$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.dec_1 = SpatialPoints(cbind(prcp_1$LON, prcp_1$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.dec_2 = SpatialPoints(cbind(prcp_2$LON, prcp_2$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))

cord.UTM_0 <- as.data.frame(spTransform(cord.dec_0, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs ")))
cord.UTM_1 <- as.data.frame(spTransform(cord.dec_1, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs ")))
cord.UTM_2 <- as.data.frame(spTransform(cord.dec_2, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs ")))

prcp_0$x1 = cord.UTM_0$coords.x1
prcp_0$x2 = cord.UTM_0$coords.x2

prcp_1$x1 = cord.UTM_1$coords.x1
prcp_1$x2 = cord.UTM_1$coords.x2

prcp_2$x1 = cord.UTM_2$coords.x1
prcp_2$x2 = cord.UTM_2$coords.x2

coordinates(prcp_0) = c('x1','x2')
coordinates(prcp_1) = c('x1','x2')
coordinates(prcp_2) = c('x1','x2')

prcp_0$dist<-sqrt((prcp_0$x1-sea$x)^2+(prcp_0$x2-sea$y)^2)
prcp_1$dist<-sqrt((prcp_1$x1-sea$x)^2+(prcp_1$x2-sea$y)^2)
prcp_2$dist<-sqrt((prcp_2$x1-sea$x)^2+(prcp_2$x2-sea$y)^2)


hist(prcp_0$value, breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')
hist(prcp_1$value, breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')
hist(prcp_2$value, breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')

hist(log(prcp_0$value), breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')
hist(log(prcp_1$value), breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')
hist(log(prcp_2$value), breaks=16, col="grey", main='Histogram of Zn', prob = TRUE, xlab = 'Zn')



xyplot(log(value) ~ dist, as.data.frame(prcp_0))
xyplot(log(value) ~ log(ELEV), as.data.frame(prcp_0))
xyplot(value ~ ELEV + dist, as.data.frame(prcp_0))

xyplot(value ~ dist, as.data.frame(prcp_1))
xyplot(value ~ ELEV, as.data.frame(prcp_1))
xyplot(value ~ ELEV + dist, as.data.frame(prcp_1))

xyplot(value ~ dist, as.data.frame(prcp_2))
xyplot(value ~ ELEV, as.data.frame(prcp_2))
xyplot(value ~ ELEV + dist, as.data.frame(prcp_2))


fit_list <- c()
#############
# no effect

v.dir_0 <- variogram(value~1,prcp_0,alpha=(0:7)*22)
v.anis_0 <- vgm(15, "Gau", 550, 5, anis=c(88, 0.9))
plot(v.dir_0, v.anis_0, pch=19)

v.fit_0 = fit.variogram(v.dir_0, v.anis_0)
fit_list[[1]] <- v.fit_0

attr(v.fit_0, 'SSErr')  # 390.1311

#############
# no effect

v.dir_1 <- variogram(value~1,prcp_1,alpha=(0:7)*22)
v.anis_1 <- vgm(15, "Gau", 550, 5, anis=c(88, 0.6))
plot(v.dir_1, v.anis_1, pch=19)

v.fit_1 = fit.variogram(v.dir_1, v.anis_1)
fit_list[[2]] <- v.fit_1

attr(v.fit_1, 'SSErr')  #  619.958


#############
# no effect

v.dir_2 <- variogram(value~1,prcp_2,alpha=(0:7)*22)
v.anis_2 <- vgm(40, "Gau", 550, 12, anis=c(88, 0.8))
plot(v.dir_2, v.anis_2, pch=19)

v.fit_2 = fit.variogram(v.dir_2, v.anis_2)
fit_list[[3]] <- v.fit_2

attr(v.fit_2, 'SSErr')  # 1747.125



###########################################################

krig.cv_0 <- krige.cv(value ~ 1, prcp_0, v.fit_0)
krig.cv_1 <- krige.cv(value ~ 1, prcp_1, v.fit_1)
krig.cv_2 <- krige.cv(value ~ 1, prcp_2, v.fit_2)

krig.cv_mse_0 <- round(sqrt(mean(krig.cv_0$residual^2)), 2)   # 2.21
krig.cv_mse_1 <- round(sqrt(mean(krig.cv_1$residual^2)), 2)   #  1.8
krig.cv_mse_2 <- round(sqrt(mean(krig.cv_2$residual^2)), 2)   # 3.65


###########################################################







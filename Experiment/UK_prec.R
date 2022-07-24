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



xyplot(value ~ sqrt(dist), as.data.frame(prcp_0))
xyplot(value ~ sqrt(ELEV), as.data.frame(prcp_0))
xyplot(value ~ log(ELEV), as.data.frame(prcp_0))
xyplot(value ~ ELEV + dist, as.data.frame(prcp_0))
# maybe dist, sì sqrt(elev), maybe log(ELEV)

xyplot(value ~ dist, as.data.frame(prcp_1))
xyplot(value ~ ELEV, as.data.frame(prcp_1))
xyplot(value ~ sqrt(ELEV), as.data.frame(prcp_1))
xyplot(value ~ log(ELEV), as.data.frame(prcp_1))
xyplot(value ~ ELEV + dist, as.data.frame(prcp_1))
# sqrt(ELEV), no dist

xyplot(value ~ dist, as.data.frame(prcp_2))
xyplot(value ~ sqrt(dist), as.data.frame(prcp_2))
xyplot(value ~ log(dist), as.data.frame(prcp_2))
xyplot(value ~ sqrt(ELEV), as.data.frame(prcp_2))
xyplot(value ~ log(ELEV), as.data.frame(prcp_2))
xyplot(value ~ ELEV, as.data.frame(prcp_2))
xyplot(value ~ ELEV + dist, as.data.frame(prcp_2))
# sqrt(ELEV)

fit_list <- c()
#############
# ELEV

v.dir_0 <- variogram(value~ELEV,prcp_0,alpha=(0:7)*22)
v.anis_0 <- vgm(22, "Gau", 550, 5, anis=c(0, 0.9))
plot(v.dir_0, v.anis_0, pch=19)

v.fit_0 = fit.variogram(v.dir_0, v.anis_0)
fit_list[[1]] <- v.fit_0

attr(v.fit_0, 'SSErr')  # 30.31595

#############
# sqrt(ELEV)

v.dir_0b <- variogram(value~sqrt(ELEV),prcp_0,alpha=(0:7)*22)
v.anis_0b <- vgm(10, "Gau", 550, 4, anis=c(132, 0.9))
plot(v.dir_0b, v.anis_0b, pch=19)

v.fit_0b = fit.variogram(v.dir_0b, v.anis_0b)
fit_list[[2]] <- v.fit_0b

attr(v.fit_0b, 'SSErr')  # 23.02392

#############
# log(ELEV)

v.dir_0c <- variogram(value~log(ELEV),prcp_0,alpha=(0:7)*22)
v.anis_0c <- vgm(22, "Gau", 550, 5, anis=c(132, 0.9))
plot(v.dir_0c, v.anis_0c, pch=19)

v.fit_0c = fit.variogram(v.dir_0c, v.anis_0c)
fit_list[[3]] <- v.fit_0c

attr(v.fit_0c, 'SSErr')  # 28.54018

#############
# sqrt(dist)

v.dir_0d <- variogram(value~sqrt(dist),prcp_0,alpha=(0:7)*22)
v.anis_0d <- vgm(10, "Gau", 550, 5, anis=c(154, 0.9))
plot(v.dir_0d, v.anis_0d, pch=19)

v.fit_0d = fit.variogram(v.dir_0d, v.anis_0d)
fit_list[[4]] <- v.fit_0d

attr(v.fit_0d, 'SSErr')  #  299.600

#############
# dist

v.dir_0e <- variogram(value~dist,prcp_0,alpha=(0:7)*22)
v.anis_0e <- vgm(10, "Gau", 550, 5, anis=c(154, 0.9))
plot(v.dir_0e, v.anis_0e, pch=19)

v.fit_0e = fit.variogram(v.dir_0e, v.anis_0e)
fit_list[[5]] <- v.fit_0e

attr(v.fit_0e, 'SSErr')  #  278.549

#############
# dist + ELEV

v.dir_0f <- variogram(value~ELEV + dist,prcp_0,alpha=(0:7)*22)
v.anis_0f <- vgm(20, "Gau", 550, 5, anis=c(0, 0.9))
plot(v.dir_0f, v.anis_0f, pch=19)

v.fit_0f = fit.variogram(v.dir_0f, v.anis_0f)
fit_list[[6]] <- v.fit_0f

attr(v.fit_0f, 'SSErr')  #  27.80117



#############
# ELEV

v.dir_1 <- variogram(value~ELEV,prcp_1,alpha=(0:7)*22)
v.anis_1 <- vgm(22, "Gau", 550, 4, anis=c(154, 0.9))
plot(v.dir_1, v.anis_1, pch=19)

v.fit_1 = fit.variogram(v.dir_1, v.anis_1)
fit_list[[7]] <- v.fit_1

attr(v.fit_1, 'SSErr')  #  154.6512

#############
# sqrt(ELEV)

v.dir_1b <- variogram(value~sqrt(ELEV),prcp_1,alpha=(0:7)*22)
v.anis_1b <- vgm(13, "Gau", 550, 5, anis=c(132, 0.9))
plot(v.dir_1b, v.anis_1b, pch=19)

v.fit_1b = fit.variogram(v.dir_1b, v.anis_1b)
fit_list[[8]] <- v.fit_1b

attr(v.fit_1b, 'SSErr')  #  157.9072

#############
# log(ELEV)

v.dir_1c <- variogram(value~log(ELEV),prcp_1,alpha=(0:7)*22)
v.anis_1c <- vgm(25, "Gau", 550, 8, anis=c(132,0.9))
plot(v.dir_1c, v.anis_1c, pch=19)

v.fit_1c = fit.variogram(v.dir_1c, v.anis_1c)
fit_list[[9]] <- v.fit_1c

attr(v.fit_1c, 'SSErr')  #  172.878

#############
# ELEV + dist

v.dir_1d <- variogram(value~ELEV + dist,prcp_1,alpha=(0:7)*22)
v.anis_1d <- vgm(30, "Gau", 550, 8, anis=c(132,0.9))
plot(v.dir_1d, v.anis_1d, pch=19)

v.fit_1d = fit.variogram(v.dir_1d, v.anis_1d)
fit_list[[10]] <- v.fit_1d

attr(v.fit_1d, 'SSErr')  #  77.47647

#############
# dist

v.dir_1e <- variogram(value~dist + dist,prcp_1,alpha=(0:7)*22)
v.anis_1e <- vgm(17, "Gau", 550, 5, anis=c(154,0.9))
plot(v.dir_1e, v.anis_1e, pch=19)

v.fit_1e = fit.variogram(v.dir_1e, v.anis_1e)
fit_list[[13]] <- v.fit_1e

attr(v.fit_1e, 'SSErr')  #  700.0778

#############
# sqrt(dist)

v.dir_1f <- variogram(value~sqrt(dist) + dist,prcp_1,alpha=(0:7)*22)
v.anis_1f <- vgm(25, "Gau", 550, 5, anis=c(132,0.9))
plot(v.dir_1f, v.anis_1f, pch=19)

v.fit_1f = fit.variogram(v.dir_1f, v.anis_1f)
fit_list[[14]] <- v.fit_1f

attr(v.fit_1f, 'SSErr')  #  504.4037




#############
# sqrt(ELEV)

v.dir_2 <- variogram(value~sqrt(ELEV),prcp_2,alpha=(0:7)*22)
v.anis_2 <- vgm(20, "Gau", 550, 14, anis=c(110, 0.9))
plot(v.dir_2, v.anis_2, pch=19)

v.fit_2 = fit.variogram(v.dir_2, v.anis_2)
fit_list[[11]] <- v.fit_2

attr(v.fit_2, 'SSErr')  # 698.973

#############
# log(ELEV)

v.dir_2b <- variogram(value~log(ELEV),prcp_2,alpha=(0:7)*22)
v.anis_2b <- vgm(35, "Gau", 550, 16, anis=c(132, 0.9))
plot(v.dir_2b, v.anis_2b, pch=19)

v.fit_2b = fit.variogram(v.dir_2b, v.anis_2b)
fit_list[[12]] <- v.fit_2b

attr(v.fit_2b, 'SSErr')  # 772.9253

#############
# sqrt(dist)

v.dir_2c <- variogram(value~sqrt(dist),prcp_2,alpha=(0:7)*22)
v.anis_2c <- vgm(25, "Gau", 550, 12, anis=c(132, 0.9))
plot(v.dir_2c, v.anis_2c, pch=19)

v.fit_2c = fit.variogram(v.dir_2c, v.anis_2c)
fit_list[[15]] <- v.fit_2c

attr(v.fit_2c, 'SSErr')  # 1881.106



###########################################################

krig.cv_0 <- krige.cv(value ~ ELEV, prcp_0, v.fit_0)
krig.cv_0b <- krige.cv(value ~ sqrt(ELEV), prcp_0, v.fit_0b)
krig.cv_0c <- krige.cv(value ~ log(ELEV), prcp_0, v.fit_0c)
krig.cv_0d <- krige.cv(value ~ sqrt(dist), prcp_0, v.fit_0d)
krig.cv_0e <- krige.cv(value ~ dist, prcp_0, v.fit_0e)
krig.cv_0f <- krige.cv(value ~ ELEV + dist, prcp_0, v.fit_0f)

krig.cv_1 <- krige.cv(value ~ ELEV, prcp_1, v.fit_1)
krig.cv_1b <- krige.cv(value ~ sqrt(ELEV), prcp_1, v.fit_1b)
krig.cv_1c <- krige.cv(value ~ log(ELEV), prcp_1, v.fit_1c)
krig.cv_1d <- krige.cv(value ~ ELEV + dist, prcp_1, v.fit_1d)
krig.cv_1e <- krige.cv(value ~ dist, prcp_1, v.fit_1e)
krig.cv_1f <- krige.cv(value ~ sqrt(dist), prcp_1, v.fit_1f)

krig.cv_2 <- krige.cv(value ~ sqrt(ELEV), prcp_2, v.fit_2)
krig.cv_2b <- krige.cv(value ~ log(ELEV), prcp_2, v.fit_2b)
krig.cv_2c <- krige.cv(value ~ sqrt(dist), prcp_2, v.fit_2c)

krig.cv_mse_0 <- round(sqrt(mean(krig.cv_0$residual^2)), 2)     # 1.57
krig.cv_mse_0b <- round(sqrt(mean(krig.cv_0b$residual^2)), 2)   # 1.61
krig.cv_mse_0c <- round(sqrt(mean(krig.cv_0c$residual^2)), 2)   # 1.72
krig.cv_mse_0d <- round(sqrt(mean(krig.cv_0d$residual^2)), 2)   # 1.66
krig.cv_mse_0e <- round(sqrt(mean(krig.cv_0e$residual^2)), 2)   # 1.64
krig.cv_mse_0f <- round(sqrt(mean(krig.cv_0f$residual^2)), 2)   # 1.55

krig.cv_mse_1 <- round(sqrt(mean(krig.cv_1$residual^2)), 2)       # 2.02
krig.cv_mse_1b <- round(sqrt(mean(krig.cv_1b$residual^2)), 2)     # 2.05
krig.cv_mse_1c <- round(sqrt(mean(krig.cv_1c$residual^2)), 2)     # 2.14
krig.cv_mse_1d <- round(sqrt(mean(krig.cv_1d$residual^2)), 2)     # 2
krig.cv_mse_1e <- round(sqrt(mean(krig.cv_1e$residual^2)), 2)     # 2.11
krig.cv_mse_1f <- round(sqrt(mean(krig.cv_1f$residual^2)), 2)     # 2.09

krig.cv_mse_2 <- round(sqrt(mean(krig.cv_2$residual^2)), 2)       # 3.57
krig.cv_mse_2b <- round(sqrt(mean(krig.cv_2b$residual^2)), 2)     # 3.6
krig.cv_mse_2c <- round(sqrt(mean(krig.cv_2c$residual^2)), 2)

###########################################################

#  resume MSE
#
#               67-73    87-93     07_13
#
#          OK:   1.8     2.21       3.65
#   UK:  ELEV:   *1.57   *2.02      *3.57
#  sqrt(ELEV):   1.61    2.05       3.61
#   log(ELEV):   1.72    2.14       3.6
#        dist:   1.64    2.11       3.66
#  sqrt(dist):   1.66    2.09       3.63
# dist + ELEV:   1.55      2        3.6       





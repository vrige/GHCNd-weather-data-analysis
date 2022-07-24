setwd("C:/Users/enogj/Desktop/proj stat/scripts/Kriging")

library(data.table)
library(dplyr) 
library(tidyr)

library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics
library(rgdal)
library(sf)

library(ggplot2)
library(maps)
library(mapproj)

## Functions for graphics 
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}

setwd("C:/Users/enogj/Desktop/proj stat/scripts/Kriging/RandomFields/R")
for (i in list.files()){source(i)}
setwd("C:/Users/enogj/Desktop/proj stat/scripts/Kriging/geoR/R")
for (i in list.files()){source(i)}

# Units of measure are tenths of degrees C
tmax = TMAX_67.73_mean_krg
tmax$value = tmax$value/10

tmax_1 = TMAX_87.93_mean_krg
tmax_1$value = tmax_1$value/10

tmax_2 = TMAX_07.13_mean_krg
tmax_2$value = tmax_2$value/10

# Convert LAT and LON to UTM
cord.dec = SpatialPoints(cbind(tmax$LON, tmax$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.dec1 = SpatialPoints(cbind(tmax_1$LON, tmax_1$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
cord.dec2 = SpatialPoints(cbind(tmax_2$LON, tmax_2$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))

#cord.UTM <- spTransform(cord.dec, CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs "))

cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
cord.UTM1 <- spTransform(cord.dec1, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
cord.UTM2 <- spTransform(cord.dec2, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))

cord.UTM = as.data.frame(cord.UTM)
cord.UTM1 = as.data.frame(cord.UTM1)
cord.UTM2 = as.data.frame(cord.UTM2)

tmax$x1 = cord.UTM$coords.x1
tmax$x2 = cord.UTM$coords.x2

tmax_1$x1 = cord.UTM1$coords.x1
tmax_1$x2 = cord.UTM1$coords.x2

tmax_2$x1 = cord.UTM2$coords.x1
tmax_2$x2 = cord.UTM2$coords.x2

# Point for distance from sea
sea=data.frame(LAT=27.5, LON=-96) 
cord.dec = SpatialPoints(cbind(sea$LON, sea$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))

cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
cord.UTM = as.data.frame(cord.UTM)

sea$x = cord.UTM$coords.x1
sea$y = cord.UTM$coords.x2

coordinates(sea)=c('x','y')

# Create a spatial dataset

coordinates(tmax) = c('x1','x2')
coordinates(tmax_1) = c('x1','x2')
coordinates(tmax_2) = c('x1','x2')

# Adding distance from sea
tmax$dist<-sqrt((tmax$x1-sea$x)^2+(tmax$x2-sea$y)^2)
tmax_1$dist<-sqrt((tmax_1$x1-sea$x)^2+(tmax_1$x2-sea$y)^2)
tmax_2$dist<-sqrt((tmax_2$x1-sea$x)^2+(tmax_2$x2-sea$y)^2)

# Exploratory analysis

hist(tmax$value, breaks=16, col="grey", main='Histogram of Tmax 67-73', prob = TRUE, xlab = 'Zn')
hist(tmax_1$value, breaks=16, col="grey", main='Histogram of Tmax 87-93', prob = TRUE, xlab = 'Zn')
hist(tmax_2$value, breaks=16, col="grey", main='Histogram of Tmax 2007-2013', prob = TRUE, xlab = 'Zn')

xyplot(value ~ dist, as.data.frame(tmax))
xyplot(value ~ log(ELEV), as.data.frame(tmax))
xyplot(value ~ ELEV + dist, as.data.frame(tmax))

xyplot(value ~ dist, as.data.frame(tmax_1))
xyplot(value ~ log(ELEV), as.data.frame(tmax_1))
xyplot(value ~ ELEV + dist, as.data.frame(tmax_1))

xyplot(value ~ sqrt(dist), as.data.frame(tmax_2))
xyplot(value ~ log(ELEV), as.data.frame(tmax_2))
xyplot(value ~ ELEV + dist, as.data.frame(tmax_2))


# Estimate the variogram

svgm <- variogram(value ~ 1, tmax)
plot(svgm, main = 'Sample Variogram tmax 65-75',pch=19,)

svgm <- variogram(value ~ ELEV, tmax)
plot(svgm, main = 'Sample Variogram tmax 05-15',pch=19,)

svgm <- variogram(value ~ ELEV+dist, tmax)
plot(svgm, main = 'Sample Variogram tmax 05-15',pch=19,)

# Look for anisotropy

plot(variogram(value ~ 1, tmax, alpha = c(0, 45, 90, 135)),main = "Var TMAX 07-13 1",pch=19)
plot(variogram(value ~ ELEV, tmax, alpha = c(0, 45, 90, 135)),main = "Var TMAX 07-13 elev",pch=19)
plot(variogram(value ~ ELEV+dist, tmax, alpha = c(0, 45, 90, 135)),main = "Var TMAX 07-13 elev+dist",pch=19)

# Select a model

vgm() # List of models

v.dir <- variogram(value~1,tmax,alpha=(0:3)*45)
v.anis <- vgm(4, "Gau", 550, 0.5, anis=c(90, 0.6))

print(plot(v.dir, v.anis, pch=19))

v.diruk <- variogram(value~ELEV,tmax,alpha=(0:3)*45)
v.anisuk <- vgm(3, "Gau", 650, 0.34, anis=c(135, 0.6))

print(plot(v.diruk, v.anisuk, pch=19))

v.dir1 <- variogram(value~1,tmax_1,alpha=(0:3)*45)
v.anis1 <- vgm(4, "Gau", 550, 0.5, anis=c(135, 0.5))

print(plot(v.dir1, v.anis1, pch=19))

v.dir1uk <- variogram(value~ELEV,tmax_1,alpha=(0:3)*45)
v.anis1uk <- vgm(3, "Gau", 600, 0.6, anis=c(135, 0.65))

print(plot(v.dir1uk, v.anis1uk, pch=19))

v.dir2 <- variogram(value~1,tmax_2,alpha=(0:3)*45)
v.anis2 <- vgm(3.5, "Gau", 500, 0.5, anis=c(90, 0.6))

print(plot(v.dir2, v.anis2, pch=19))

v.dir2uk <- variogram(value~ELEV,tmax_2,alpha=(0:3)*45)
v.anis2uk <- vgm(3, "Gau", 600, 0.5, anis=c(135, 0.65))

print(plot(v.dir2uk, v.anis2uk, pch=19))

# Fit the model
# test = fit.variogram.reml(value~dist+ELEV, tmax, model=v.anis)
# test

v.fit = fit.variogram(v.dir, model=v.anis)
v.fit

v.fituk = fit.variogram(v.diruk, model=v.anisuk)
v.fituk

v.fit1 = fit.variogram(v.dir1, model=v.anis1)
v.fit1

v.fit1uk = fit.variogram(v.dir1uk, model=v.anis1uk,fit.ranges = TRUE)
v.fit1uk

v.fit2 = fit.variogram(v.dir2, model=v.anis2,fit.sills = FALSE,fit.ranges = TRUE)
v.fit2

v.fit2uk = fit.variogram(v.dir2uk, model=v.anis2uk,fit.sills = TRUE,fit.ranges = FALSE)
v.fit2uk
# test = fit.variogram.gls(value~dist+ELEV,data=tmax, model=v.anis,plot=TRUE)

print(plot(v.dir2uk,v.fit2uk,pch=19))


print(plot(svgm, test, pch=19))

#############################################################
##############                                 ##############
######          SPATIAL PREDICTION & KRIGING          #######
##############                                 ##############
#############################################################

hist(tmax$value, breaks=16, col="grey", main='Histogram of Tmax', prob = TRUE, xlab = 'Tmax')


# Visualization
library(ggplot2)
library(maps)
library(mapproj)

head(map_data('state'),15)

usa_tbl <- map_data("state", region = c('texas')) %>% as_tibble()

# usa_tbl %>% 
#   ggplot() + 
#   geom_map(
#     map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
#   ) +
#   coord_map("ortho",orientation = c(39,-98,0))+
#   geom_point(data = TMAX_07.13_mean_krg, aes(x=LON,y=LAT,color = value), size = 3) +
#   scale_color_gradient(low="blue", high="yellow",name="TMAX")

# Kriging 

# s0.new=data.frame(LAT=33.67407, LON=-100.43701, ELEV=601) 
# cord.dec = SpatialPoints(cbind(s0.new$LON, s0.new$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
# 
# cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
# cord.UTM = as.data.frame(cord.UTM)
# 
# s0.new$x = cord.UTM$coords.x1
# s0.new$y = cord.UTM$coords.x2
# 
# 
# coordinates(s0.new)=c('x','y','ELEV')
# 
# g.tr <- gstat(formula = value ~ 1, data = tmax, model = test)
# 
# predict(g.tr, s0.new, BLUE = FALSE)



# CREATING A GRID FOR PREDICTION
# Imported dataset convert to spatial
grid

colnames(grid) = c('x','y','ELEV','LAT','LON')

coordinates(grid)=c('x','y')

# Adding distance from sea

grid$dist<-sqrt((grid$x-sea$x)^2+(grid$y-sea$y)^2)

# 67-73

g.tr <- gstat(formula = value ~ ELEV, data = tmax, model = v.fituk)

res = predict(g.tr, grid, BLUE = FALSE)

fin = as.data.frame(res)

cord.UTM = SpatialPoints(cbind(fin$x, fin$y), proj4string = CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))

cord.dec <- spTransform(cord.UTM,CRS("+proj=longlat +datum=WGS84"))

cord.dec = as.data.frame(cord.dec)

fin$LON = cord.dec$coords.x1
fin$LAT = cord.dec$coords.x2

# 87-93

g.tr1 <- gstat(formula = value ~ ELEV, data = tmax_1, model = v.fit1uk)

res1 = predict(g.tr1, grid, BLUE = FALSE)

fin1 = as.data.frame(res1)

cord.UTM = SpatialPoints(cbind(fin1$x, fin1$y), proj4string = CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))

cord.dec <- spTransform(cord.UTM,CRS("+proj=longlat +datum=WGS84"))

cord.dec = as.data.frame(cord.dec)

fin1$LON = cord.dec$coords.x1
fin1$LAT = cord.dec$coords.x2

# 2007-2013

g.tr2 <- gstat(formula = value ~ ELEV, data = tmax_2, model = v.fit2uk)

res2 = predict(g.tr2, grid, BLUE = FALSE)

fin2 = as.data.frame(res2)

cord.UTM = SpatialPoints(cbind(fin2$x, fin2$y), proj4string = CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))

cord.dec <- spTransform(cord.UTM,CRS("+proj=longlat +datum=WGS84"))

cord.dec = as.data.frame(cord.dec)

fin2$LON = cord.dec$coords.x1
fin2$LAT = cord.dec$coords.x2

# CROSS VALIDATION

krig.cv <- krige.cv(value ~ 1, tmax, v.fit)
krig.cvuk <- krige.cv(value ~ sqrt(dist), tmax, v.fituk)

krig.cv1 <- krige.cv(value ~ 1, tmax_1, v.fit1)
krig.cv1uk <- krige.cv(value ~ sqrt(dist),  tmax_1, v.fit1uk)

krig.cv2 <- krige.cv(value ~ 1,  tmax_2, v.fit2)
krig.cv2uk <- krige.cv(value ~ sqrt(dist),  tmax_2, v.fit2uk)


krig.cv_mse <- round(sqrt(mean(krig.cv$residual^2)), 2)   # 0.89
krig.cv_mseuk <- round(sqrt(mean(krig.cvuk$residual^2)), 2)   # 0.85

krig.cv_mse1 <- round(sqrt(mean(krig.cv1$residual^2)), 2)   # 1.05
krig.cv_mse1uk <- round(sqrt(mean(krig.cv1uk$residual^2)), 2)   # 1.04

krig.cv_mse2 <- round(sqrt(mean(krig.cv2$residual^2)), 2)   # 1.01
krig.cv_mse2uk <- round(sqrt(mean(krig.cv2uk$residual^2)), 2)   # 1.03

# resume MSE  TMAX
#
#              67-73    87-93     07_13
#
#         OK:   0.89     1.05      1.01
#  UK:  ELEV:   0.65     0.69      0.60     
# sqrt(ELEV):   0.75     0.81      0.76
#  log(ELEV):   0.86     1.03      1.02
#       dist:   0.85     1.03      1.03
# sqrt(dist):   0.85     1.04      1.03
#  ELEV+dist:   0.64     0.69      0.63


# Plotting
head(map_data('state'),15)

usa_tbl <- map_data("state", region = c('texas')) %>% as_tibble()

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = TMAX_07.13_mean_krg, aes(x=LON,y=LAT,color = value/10), size = 3) + scale_size(name="TMAX")+
  scale_color_gradientn(colors = rev(seis), name = "°C") +
  ggtitle("Mean Maximum Temperature between 2007-2013")

# Plotting results

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin, aes(x=LON,y=LAT,color = var1.pred),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= rev(seis) ,name="°C",limits=c(18,31.25534))+
  #ggtitle("Maximum Temperature Prediction 1967-1973")+
  labs(y = "", x = "")+
  theme(plot.title = element_text(face="bold",size = 15))
  

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin1, aes(x=LON,y=LAT,color = var1.pred),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= rev(seis) ,name="°C",limits=c(18,31.25534))+
  #ggtitle("Maximum Temperature Prediction 1987-1993")+
  labs(y = "", x = "")

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin2, aes(x=LON,y=LAT,color = var1.pred),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= rev(seis) ,name="°C",limits=c(18,31.25534))+
  #ggtitle("Maximum Temperature Prediction 2007-2013")+
  labs(y = "", x = "")




usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = diff, aes(x=LON,y=LAT,color = var1.pred),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= rev(jet) ,name="°C", values = c(1,0.65,0))+
  #ggtitle("Maximum Temperature difference 70 10")+
  labs(y = "", x = "")

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = diff3, aes(x=LON,y=LAT,color = var1.pred),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= jet ,name="°C")+
  ggtitle("Maximum Temperature difference 70 90")+
  labs(y = "Latitude", x = "Longitude")


jet = c("#00007F", 
        "#0000B2", "#0000E5", "#0019FF", "#004DFF", "#007FFF", "#00B2FF", 
        "#00E5FF", "#FFFFF2", "#FFFFD9", "#FFFFBF", "#FFFFA5", "#FFFF8C", 
        "#FFE500", "#FFB300", "#FF7F00", "#FF4C00", "#FF1900", "#E50000", 
        "#B20000")






usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin, aes(x=LON,y=LAT,color = var1.var),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= rev(pal_var) ,name="Variance")+
  ggtitle("Maximum Temperature Variance 1967-1973")+
  labs(y = "Latitude", x = "Longitude")

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin1, aes(x=LON,y=LAT,color = var1.var),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= pal_var ,name="Variance",values = c(1.0,0.8,0.1,0.05,0))+
  #ggtitle("Maximum Temperature Variance 1987-1993")+
  labs(y = "", x = "")

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin2, aes(x=LON,y=LAT,color = var1.var),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= pal_var ,name="Variance",values = c(1.0,0.8,0.5,0.1,0.05,0))+
  ggtitle("Maximum Temperature Variance 2007-2013")+
  labs(y = "Latitude", x = "Longitude")





usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = fin, aes(x=LON,y=LAT,color = var1.var),size=1.8) +
  #scale_colour_stepsn(n.breaks = 20, colours = seis)
  scale_colour_gradientn(colours= rev(pal_var) ,name="Variance")+
  ggtitle("Maximum Temperature Variance 1965-1973")+
  labs(y = "Latitude", x = "Longitude")

pal_var =c("red","yellow", "aquamarine","light blue","blue")
seis = c("#AA0000", 
         "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", 
         "#FFB700", "#FFDD00", "#FFFF00", "#FFFF00", "#FFFF00", "#BDFF0C", 
         "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", 
         "#0024E3")

prj_dd <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "


df_elev_epqs <- get_elev_point(grid, prj = prj_dd, src = "epqs")

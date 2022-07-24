library(data.table)
library(dplyr) 
library(tidyr)

library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)        ## Geostatistics
library(rgdal)
library(sf)

## Functions for graphics 
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est<-function(x,C0, ...){C0-cov.spatial(x, ...)}

setwd("C:/Users/enogj/Desktop/proj stat/scripts/Kriging/RandomFields/R")
for (i in list.files()){source(i)}
setwd("C:/Users/enogj/Desktop/proj stat/scripts/Kriging/geoR/R")
for (i in list.files()){source(i)}

# Units of measure are tenths of degrees C
tmin_6773 = TMIN_67.73_mean_krg
tmin_6773$value = tmin_6773$value/10

# tmin.1 = TMIN_67.73_mean_krg
# tmin.1$value = tmin$value/10

# tmin_1 = tmin_87.93_mean_krg
# tmin_1$value = tmin_1$value/10
# 
# tmin_2 = tmin_07.13_mean_krg
# tmin_2$value = tmin_2$value/10

sea=data.frame(LAT=27.5, LON=-96) 
cord.dec = SpatialPoints(cbind(sea$LON, sea$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))

cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
cord.UTM = as.data.frame(cord.UTM)

sea$x = cord.UTM$coords.x1
sea$y = cord.UTM$coords.x2


coordinates(sea)=c('x','y')

tmin_6773$dist<-sqrt((tmin_6773$x1-sea$x)^2+(tmin_6773$x2-sea$y)^2)

# Convert LAT and LON to UTM
cord.dec = SpatialPoints(cbind(tmin_6773$LON, tmin_6773$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))

#cord.UTM <- spTransform(cord.dec, CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs "))

cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))


cord.UTM = as.data.frame(cord.UTM)

tmin_6773$x1 = cord.UTM$coords.x1
tmin_6773$x2 = cord.UTM$coords.x2

# tmin_1$x1 = cord.UTM$coords.x1
# tmin_1$x2 = cord.UTM$coords.x2

# Create a spatial dataset

coordinates(tmin_6773) = c('x1','x2')
tmin_6773$dist<-sqrt((tmin_6773$x1-sea$x)^2+(tmin_6773$x2-sea$y)^2)

# coordinates(tmin.1)=c('LON','LAT')
# coordinates(tmin_1) = c('x1','x2')

# Estimate the variogram

svgm <- variogram(value ~ dist+ELEV, tmin_6773)
plot(svgm, main = 'Sample Variogram tmin 67-73',pch=19)

# svgm <- variogram(value ~ ELEV, tmin_1)
# plot(svgm, main = 'Sample Variogram tmin 05-15',pch=19,)

# Look for anisotropy

#plot(variogram(value ~ dist+ELEV, tmin, alpha = c(0, 22, 45, 67, 90, 112, 135, 158)),main = "Variograms tmin 1967-1973",pch=19)
# plot(variogram(value ~ ELEV, tmin_1, alpha = c(0, 45, 90, 135)),main = "Variograms tmin 2015",pch=19)

# Select a model

vgm() # List of models

# v1 <- vgm(1.2, "Gau", 500, 0.35, anis=c(90, 0.5))
# v2 <- vgm(0.6,"Gau",500, 0.35, anis=c(90, 0.5), add.to = v1)
# v2


v.dir <- variogram(value~dist+ELEV,tmin_6773,alpha=c(0:3)*45)
plot(v.dir,pch=19)
v.anis <- vgm(1.5, "Gau", 500, 0.5, anis=c(90, 0.6))
#v.anis2<- vgm(0.9313724,290.6042,0.4227548,anis=c(45, 0.6))
#plot(variogram(value~dist+ELEV,tmin,alpha=c(0,90)),pch=19)

print(plot(v.dir, v.anis, pch=19))
#print(plot(v.dir, v.anis2, pch=19))

# Fit the model
test = fit.variogram.reml(value~dist+ELEV, tmin, model=v.anis)
test



test = fit.variogram(v.dir, model=v.anis)

print(plot(test,700))


print(plot(svgm, test, pch=19))
test

#############################################################
##############                                 ##############
######          SPATIAL PREDICTION & KRIGING          #######
##############                                 ##############
#############################################################

hist(tmin$value, breaks=16, col="grey", main='Histogram of tmin', prob = TRUE, xlab = 'tmin')


# Visualization
library(ggplot2)
library(maps)
library(mapproj)

head(map_data('state'),15)

usa_tbl <- map_data("state", region = c('texas')) %>% as_tibble()

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = TMIN_67.73_mean_krg, aes(x=LON,y=LAT,color = value), size = 5) + scale_size(name="tmin")+
  scale_color_gradient(low="blue", high="yellow")

# Kriging 

# s0.new=data.frame(LAT=33.67407, LON=-100.43701, ELEV=601) 
# cord.dec = SpatialPoints(cbind(s0.new$LON, s0.new$LAT), proj4string = CRS("+proj=longlat +datum=WGS84"))
# 
# cord.UTM <- spTransform(cord.dec, CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))
# cord.UTM = as.data.frame(cord.UTM)
# 
# s0.new$x = cord.UTM$coords.x1
# s0.new$y = cord.UTM$coords.x2
# s0.new$dist= sqrt((s0.new$x-sea$x)^2+(s0.new$y-sea$y)^2)
# 
# 
# coordinates(s0.new)=c('x','y')

g.tr_tmin_6773 <- gstat(formula = value ~ dist+ELEV, data = tmin_6773, model = test)

# predict(g.tr, s0.new, BLUE = FALSE)
# 
# predict(g.tr, s0.new, BLUE = TRUE)



# res <- 15 # resolution
# # round extremes to resolution
# (x.min <- bbox(tmin)[1,1]%/%res*res)
# (y.min <- bbox(tmin)[2,1]%/%res*res)
# (x.max <- (bbox(tmin)[1,2]+res)%/%res*res) # make sure it is outside the bbox
# (y.max <- (bbox(tmin)[2,2]+res)%/%res*res) # make sure it is outside the bbox
# 
# x.min = 373.1241 -25
# x.max = 1592.729 +50
# 
# y.min = 423.7076 -25
# y.max = 1615.424 +25
# 
# grid <- expand.grid(x = seq(x.min, x.max, by=res),
#                     y = seq(y.min, y.max, by=res))
# grid$dist<- sqrt((grid$x-sea$x)^2+(grid$y-sea$y)^2)
# grid$ELEV<- df_elev_epqs$elevation
# 
# 
# 
# real.grid <- na.omit(grid)   # grid[-which(grid$ELEV=='NA'),]
# 
# coordinates(real.grid) <- c("x", "y")
# gridded(real.grid) <- T; fullgrid(real.grid) <- T

prova_tmin_6773 = as.data.frame(predict(g.tr_tmin_6773, grid, BLUE = FALSE))

spplot(prova)

# convert coords

cord.UTM = SpatialPoints(cbind(prova_tmin_6773$x, prova_tmin_6773$y), proj4string = CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))

cord.dec <- spTransform(cord.UTM,CRS("+proj=longlat +datum=WGS84"))

cord.dec = as.data.frame(cord.dec)

prova_tmin_6773$LON = cord.dec$coords.x1
prova_tmin_6773$LAT = cord.dec$coords.x2


# Plotting
head(map_data('state'),15)

usa_tbl <- map_data("state", region = c('texas')) %>% as_tibble()

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = prova_tmin_6773, aes(x=LON,y=LAT,color = var1.var), size = 1) + scale_size(name="tmin")+
  scale_color_gradient(low="blue", high="yellow")


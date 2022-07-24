library(raster)
library(rgdal)
library(sp)
#Use this line to download the County Shapefile for the US
us<-getData('GADM', country='USA', level=1)  
texas <- subset(us,NAME_1=="Texas")
plot(texas)

states <- map_data("state")
tx_county <- subset(states, region == "texas")

writeOGR(obj=texas, dsn="tempdir", layer="texas", driver="ESRI Shapefile")  

grid = grid[texas,]
plot(grid)


spdf <- SpatialPointsDataFrame(cbind(grid$LON,grid$LAT), data =grid,
      proj4string = CRS("+proj=longlat +datum=WGS84"))

spdf = spdf[texas,]

plot(spdf)

final_grid = as.data.frame(spdf)

plot(texas)
image(final_grid,add=T)

usa_tbl %>% 
  ggplot() + 
  geom_map(
    map = usa_tbl, aes(x=long, y=lat, map_id = region),color = "gray80", fill = "white", size = 0.4
  ) +
  coord_map("ortho",orientation = c(39,-98,0))+
  geom_point(data = final_grid, aes(x=LON,y=LAT,color = elev),size = 1.7) +
  scale_colour_gradientn(colours = c("#EA5C3A","#FAF9BA", "#77B9AC"))

########################################################################################################

# Extract elevation from data
library(elevatr)

prj_dd <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "

df_elev_epqs <- get_elev_point(grid, prj = prj_dd, src = "epqs")

elev = as.data.frame(df_elev_epqs)

grid$elev = elev$elevation

cord.UTM = SpatialPoints(cbind(grid$x, grid$y), proj4string = CRS("+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "))

cord.dec <- spTransform(cord.UTM,CRS("+proj=longlat +datum=WGS84"))
cord.dec = as.data.frame(cord.dec)
grid$LAT = cord.dec$coords.x2
grid$LON = cord.dec$coords.x1


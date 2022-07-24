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


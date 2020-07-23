##
# mapping risk predicted
### function plit polygon pixels

Features <-read.csv("Landscape and host data.csv", row.names = 1, head=T)
FeaturesnoNA<-Features[complete.cases(Features), ];str(Features)
## area pixel square meters
# a <- poly_split(shp = mapa, area_pixel = 10000)
# str(a@data,2)
library(sf)
library(ggplot2)
library(sp)
library(rgdal)
library(rgeos)

coords <- FeaturesnoNA[c("Longitude", "Latitude")]


# Letting R know that these are specifically spatial coordinates
sp1 <- SpatialPoints(coords)
plot(sp)

coordinates( FeaturesnoNA ) <- c( "Longitude", "Latitude")
proj4string( FeaturesnoNA ) <- CRS( "+proj=longlat +datum=WGS84" )
distInMeters <- 2000
pc100km <- gBuffer( FeaturesnoNA, width=100*distInMeters, byid=TRUE )
plot(pc100km)
state.ll83 <- spTransform(pc100km, CRS("+proj=longlat +ellps=32724"))
plot(state.ll83)

size <- 5000
div_5 <- poly_split(shp = state.ll83, area_pixel = size,
                    CRS_utm = "+init=epsg:32617")



ggplot(aux_df, aes(x=max(Longitude), y=max(Latitude))) + 
         geom_polygon(color='black', fill=NA) +
         coord_map() +
         theme_classic()
       
size <- 5000
div_5 <- poly_split(shp = buf, area_pixel = size,
                    CRS_utm = "+init=epsg:32617")


# dataset <- gtas_t1
# var <- "dg"
# n_removed_nodes <- 1
# total_removed <- 2000
#install.packages("progress")


pred <- raster:: predict(predictors, rf.fit,  type='prob', ext=EU)
plot(pred)

poly_split <- function(shp, area_pixel = 1000, CRS_utm = "+init=epsg:32724"){
  
  shp <- spTransform(shp, CRSobj = CRS_utm)
  
  if(length(area_pixel) == 1){
    cs <- c(1, 1)*area_pixel 
  } else {
    cs <- area_pixel
  }
  
  
  
  grdpts <- makegrid(shp, cellsize = cs, pretty = F)
  spgrd <- SpatialPoints(grdpts, proj4string = CRS(proj4string(shp)))
  spgrdWithin <- SpatialPixels(spgrd[shp,])
  spgrdWithin <- as(spgrdWithin, "SpatialPolygons")
  
  # Extract polygon ID's
  pid <- sapply(slot(spgrdWithin, "polygons"), function(x) slot(x, "ID")) 
  p.df <- data.frame( ID=1:length(spgrdWithin), ID2 = pid, row.names = pid) 
  spgrdWithin <- SpatialPolygonsDataFrame(spgrdWithin, p.df)
  spgrdWithin$ID2 <- as.character(spgrdWithin$ID2)
  
  #print(plot(spgrdWithin))
  
  return(spgrdWithin)
}

## area pixel square meters
# a <- poly_split(shp = mapa, area_pixel = 10000)
# str(a@data,2)

sfd = st_as_sf(sfd, coords=c("long","lat"), crs=27700)

size <- 5000
div_5 <- poly_split(shp = nc, area_pixel = size,
                    CRS_utm = "+init=epsg:32617")


# dataset <- gtas_t1
# var <- "dg"
# n_removed_nodes <- 1
# total_removed <- 2000
#install.packages("progress")


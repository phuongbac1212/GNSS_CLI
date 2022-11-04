require(sp)
require(sf)

location.filter <- function(shape.file = "roi", point.dataframe) {
  
  pol = read_sf(dsn = shp_path, layer = shape.file)
  pol.coor = as.data.frame(st_coordinates(pol)[,c("X", "Y") ])
  
  point = point.in.polygon(point.dataframe$Lon_simu, point.dataframe$Lat_simu, pol.coor$X, pol.coor$Y)
  index = which(point ==1)
  
  result  = point.dataframe[index,]

    # require(leaflet)
    # leaflet() %>% 
    #   addTiles() %>% 
    #   addCircleMarkers(lng = result$Lon_simu, lat = result$Lat_simu, radius = 1, color="red") %>% 
    #   addPolygons(lng= pol.coor$X, lat= pol.coor$Y)
    # 
  return(result)
}

#leaflet() %>% addTiles() %>% addCircleMarkers(lng = temp$Lon_simu, lat= temp$Lat_simu, radius = 1) %>% addPolygons(lng = pol.coor$X, lat = pol.coor$Y)



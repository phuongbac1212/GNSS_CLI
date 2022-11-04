require(rgdal)

dataframe2txt <- function(x, name){
  print("Export to 3.output/tmp/ [Plain text format]")
  write.table(x,row.names = F, file = paste0("3.output/tmp/Result-",name,".txt"))
}

dataframe2kml <- function(x, path) {
  print("Export to 3.output/ [KML format]")
  pnr = unique(res$PNR)
  lapply(
    pnr,
    FUN = function(x) {
      temp = filter(res, PNR == x) %>% drop_na()
      coordinates(temp) <-
        c("Lon_simu", "Lat_simu")          #Build a SpatialPointsData Frame
      proj4string(temp) <- CRS("+proj=longlat +datum=WGS84")
      try(writeOGR(
        temp,
        dsn = paste0(path, x, ".KML"),
        layer = x,
        driver = "KML"
      ),
      silent = T)
    }
  )
}

dataframe2map <- function(x) {
  require(leaflet)
  leaflet() %>% addTiles() %>% addCircleMarkers(lng = x$Lon_simu, lat = x$Lat_simu)
}

result.concat <- function() {
  list.result =list.files("3.output/tmp/")
  tmp = read.table(paste0("3.output/tmp/", list.result[1]),header=T) %>% drop_na()
  unlink("3.output/Result.txt")
  write.table(tmp, file = "3.output/Result.txt", append=T, col.names = T,  row.names = F)
  
  for (file in list.result[2:length(list.result)]) {
    tmp = read.table(paste0("3.output/tmp/", file),header=T) %>% drop_na() 
    write.table(tmp, file = "3.output/Result.txt", append=T, col.names = F,  row.names = F)
  }
}

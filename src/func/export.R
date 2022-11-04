require(rgdal)
require(tidyr)
dataframe2txt <- function(x, name, append = F) {
  if (nrow(x) == 0)
    return(NULL)
  if (!dir.exists(paste0(result_path, "tmp")))
    dir.create(paste0(result_path, "tmp"))
  
  file = paste0(result_path, "tmp/Result-", name, ".txt")
  print(paste("Export to", result_path,file , " [Plain text format]"))
  
  write.table(
    x,
    row.names = !file.exists(file),
    append = append,
    file = file
  )
}

dataframe2fst <- function(x, name) {
  if (nrow(x) == 0)
    return(NULL)
  print(paste(
    "Export to",
    result_path ,
    " [Fast Serialization of Data Frames format]"
  ))
  if (!dir.exists(paste0(result_path, "tmp")))
    dir.create(paste0(result_path, "tmp"))
  write_fst(x, paste0(result_path, "tmp/Result-", name, ".fst"))
}

dataframe2kml <- function(x) {
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
        dsn = paste0("3.output/", x, ".KML"),
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
  list.result = list.files(paste0(result_path, "tmp"), pattern = "txt")
  if (length(list.result) > 0) {
    tmp = read.table(paste0(result_path, "tmp/", list.result[1]), header =
                       T) %>% drop_na()
    unlink(paste0(result_path, "Result.txt"))
    write.table(
      tmp,
      file = paste0(result_path, "Result.txt"),
      append = T,
      col.names = T,
      row.names = F
    )
    
    for (file in list.result[2:length(list.result)]) {
      tmp = read.table(paste0(result_path, "tmp/", file), header = T) %>% drop_na()
      write.table(
        tmp,
        file = paste0(result_path, "Result.txt"),
        append = T,
        col.names = F,
        row.names = F
      )
    }
  }
  
  list.result.fst = list.files(paste0(result_path, "tmp"), pattern = "fst")
  if (length(list.result.fst > 0)) {
    tmp = as.data.frame(fst(paste0(
      result_path, "tmp/", list.result.fst[1]
    ))) %>% drop_na()
    if (length(list.result) == 0) {
      unlink(paste0(result_path, "Result.txt"))
    }
    write.table(
      tmp,
      file = paste0(result_path, "Result.txt"),
      append = T,
      col.names = T,
      row.names = F
    )
    
    for (file in list.result.fst[2:length(list.result.fst)]) {
      tmp = as.data.frame(fst(paste0("3.output/tmp/", file), header = T)) %>% drop_na()
      write.table(
        tmp,
        file = paste0(result_path, "Result.txt"),
        append = T,
        col.names = F,
        row.names = F
      )
    }
  }
}

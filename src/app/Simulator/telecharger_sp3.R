telecharger_sp3 <- function(fichier,path){
     
      date_cur <- strsplit(strsplit(fichier,".obs")[[1]][1], "_")[[1]]
      sp3.link  = sp3.getLink(date_cur)
      if(file.exists(sp3.getFileNamFromLink(sp3.link))) return(sp3.getFileNamFromLink(sp3.link))
      sp3.file = sp3.download(sp3.link, path)
      return(sp3.file)
}



library(RCurl)
library(lubridate)
sp3.getLink <- function(time) {
   time = as.numeric(time)
   timestamp = time[1] * 7 * 24 * 60 * 60 + (time[2] + 1) * 24 * 60 * 60
   date = as_datetime(timestamp, origin = "1980-01-06")
   week = time[1]
   year = year(date)
   day_of_year = yday(date) - 1
   temp.file = paste0(url_igs,
                      week,
                      "/WUM0MGXULA_",
                      year,
                      sprintf("%03d", day_of_year),
                      "*")
   
   print(paste("Trying to ping [", temp.file, "]"))
   
   result <-
      RCurl::getURL(
         temp.file,
         verbose = F,
         ftp.use.epsv = FALSE,
         dirlistonly = TRUE
      )
   
   result = paste(url_igs, week , "/", strsplit(result, "\r*\n")[[1]], sep = "")
   pattern.list = c(paste0("/",
                           week,
                           "/WUM0MGXULA_",
                           year,
                           sprintf("%03d", day_of_year)),
                    "_ORB.SP3.gz")
   result = result[Reduce(intersect, sapply(pattern.list, grep, result))[1]]
   print(paste("Got: [", result, "]"))
   return(result)
}

sp3.getFileNamFromLink <- function(link) {
   name = strsplit(link, "/")[[1]]
   name = strsplit(name[length(name)], ".gz")[[1]][1]
   return(name[length(name)])
}

sp3.download <- function(link, path) {
   print(paste("Start to download sp3 from ", link, "to", path))
   file.name = tail(strsplit(link, "/")[[1]], 1)
   if (file.exists(paste0(path, strsplit(file.name, ".gz")[[1]][1]))) {
      print("existed... exit...")
      return(strsplit(file.name, ".gz")[[1]][1])
   }
   download.file(link,
                 paste0(path, file.name),
                 mode = 'wb',
                 quiet = T)
   R.utils::gunzip(paste0(path, "/", file.name),
                   overwrite = T,
                   remove = T)
   return(strsplit(file.name, ".gz")[[1]])
}

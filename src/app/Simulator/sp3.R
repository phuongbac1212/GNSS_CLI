##################################################################
#
# This file contain function that download and spline sp3
# satellite product data
#
# Nov 2020 - Nguyen Phuong Bac
#
##################################################################

###############################################################

#' Read SP3 file
#'
#' @param path_input SP3 file path
#' @return A dataframe of satellite position and the corresponding time
SP3.read <- function(path_input = "sp3.fileuz") {
  sp3_p = ""
  file <- file(paste(path_input, sp3_p, sep = ''), open = "r")
  lines_p <-
    strsplit(gsub("\\s+", " ", gsub(
      "^\\s+|\\s+$", "", readLines(file, warn = FALSE)
    )), ' ')
  close(file)
  
  lines_date <-
    lines_p[substr(sapply(lines_p, '[[', 1), 1, 1) == '*']
  annee <- sapply(lines_date, '[[', 2)
  mois <- sapply(lines_date, '[[', 3)
  jour <- sapply(lines_date, '[[', 4)
  heure <- sapply(lines_date, '[[', 5)
  minute <- sapply(lines_date, '[[', 6)
  seconde <- sapply(lines_date, '[[', 7)
  
  ts <-
    mapply(
      FUN =  function(year, month, day, hour, min, sec) {
        return(as.numeric(strptime(
          paste(year, "-", month, "-", day, " ", hour, ":", min, ":", sec, sep = ""),
          format = "%Y-%m-%d %H:%M:%OS",
          tz = "+0"
        )) - 315964800)
      },
      annee,
      mois,
      jour,
      heure,
      minute,
      seconde
    )
  
  lines_sat <-
    lines_p[substr(sapply(lines_p, '[[', 1), 1, 1) == 'P']
  ID <- sapply(lines_sat, '[[', 1)
  ts <- rep(ts, each = length(unique(ID)))
  x <- sapply(lines_sat, '[[', 2)
  y <- sapply(lines_sat, '[[', 3)
  z <- sapply(lines_sat, '[[', 4)
  
  
  return(
    cbind(
      "TimeStamp" = as.numeric(ts),
      "X" = as.numeric(x) * 1000,
      "Y" = as.numeric(y) * 1000,
      "Z" = as.numeric(z) * 1000,
      "ID" = ID
    )
  )
}


require(dplyr, quietly = TRUE)
require(parallel)
source("multidimensional_spline.R")

#' Iterpolate GNSS production data
#'
#' @param path_input SP3 file path
#' @return A dataframe of satellite position and the corresponding time each 1s
SP3.spline <- function(path_input = "sp3.fileuz") {
  print("Iterpolating SP3")
  sp3.reader.df = SP3.read(path_input = path_input)
  list.sat = unique(sp3.reader.df[, "ID"])
  splined.df <- lapply(
    list.sat,
    FUN = function(sat, ...) {
      sub.df = which(sp3.reader.df[, "ID"] == sat)
      sub.df = sp3.reader.df[sub.df,]
      splined.sub.df = spline.Nd.NoPar(sub.df[, c(1, 2, 3, 4)], dimension = ncol(sub.df) - 1, step =
                                         300)
      splined.sub.df[, "ID"] <- sat
      return(splined.sub.df)
    }
  ) %>% bind_rows()
  return(splined.df)
}


#' Download LASTEST SP3 file from specified link - -OLD-
#'
#' @param url SP3 ftp link
#' @return Ftp url to the downloaded file
sp3.downloader <-
  function(url = "ftp://igs.ign.fr/pub/igs/products/mgex/") {
    unlink("sp3.fileuz")
    wn = (as.numeric(Sys.time()) - 315964782) %/% 60 %/% 60 %/% 24 %/% 7
    url = paste0(url, wn, "/")
    result <-
      RCurl::getURL(
        url,
        verbose = TRUE,
        ftp.use.epsv = FALSE,
        dirlistonly = TRUE,
        crlf = TRUE
      )
    result2 <- paste(url, strsplit(result, "\r*\n")[[1]], sep = "")
    download.file(tail(result2, n = 1),
                  "sp3.file.gz",
                  mode = 'wb',
                  quiet = TRUE)
    
    R.utils::gunzip("sp3.file.gz")
    file.rename(from = "sp3.file", to = "sp3.fileuz")
    return(tail(result2, n = 1))
  }

require(RCurl)
require(lubridate)
sp3.getLink <- function(timestamp, sp3_path = ".") {
  print("get link for sp3")
  time = as.numeric(timestamp)
  date = as_datetime(timestamp, origin = "1980-1-6")
  week = floor(time/604800)
  year = year(date)
  day_of_year = yday(date)
  temp.file = paste0(param["sp3_suffix"],
                     year,
                     sprintf("%03d", day_of_year))
  list.sp3 = list.files(path = sp3_path, pattern=param["sp3_suffix"])
  sp3.avai = list.sp3[grep(temp.file, list.sp3)]
  if(length(sp3.avai>2)) sp3.avai = sp3.avai[[1]]
  if (length(sp3.avai)>0) {
    if (length(grep(".gz", sp3.avai)) > 0) {
      R.utils::gunzip(paste0(sp3_path, sp3.avai), overwrite = T, remove = T)
    }
    return(paste0(sp3_path,sp3.avai))
  }
  temp.link = paste0(param["sp3_ftp"],
                     week,
                     "/",param["sp3_suffix"],
                     year,
                     sprintf("%03d", day_of_year),
                     "*")
  
  print(paste("Trying to ping [", temp.link, "]"))
  
  result <-
    RCurl::getURL(
      temp.link,
      verbose = F,
      ftp.use.epsv = FALSE,
      dirlistonly = TRUE
    )
  
  result = paste(param["sp3_ftp"], week , "/", strsplit(result, "\r*\n")[[1]], sep = "")
  pattern.list = c(paste0("/",
                          week,
                          "/", param["sp3_suffix"],
                          year,
                          sprintf("%03d", day_of_year)),
                   "_ORB.SP3.gz")
  result = result[Reduce(intersect, sapply(pattern.list, grep, result))[1]]
  print(paste("Got: [", result, "]"))
  return(result)
}

sp3.getFileNamFromLink <- function(link) {
  name = strsplit(link, "/")[[1]]
  return(name[length(name)])
}

sp3.download <- function(link, path) {
  print(paste("Start to download sp3 from ", link, "to", path))
  file.name = tail(strsplit(link, "/")[[1]], 1)
  # if (file.exists(paste0(sp3_path, strsplit(file.name, ".gz")[[1]][1]))) {
  #   print("existed... exit...")
  #   return(strsplit(file.name, ".gz")[[1]][1])
  # }
  xfun::download_file(link,
                paste0(path,"/", file.name),
                mode = 'wb',
                quiet = F)
  R.utils::gunzip(paste0(path, "/", file.name),
                  overwrite = T,
                  remove = T)
  return(strsplit(file.name, ".gz")[[1]])
}

# 
# #####################################
# # for download WOMONXULA from cddis
# #####################################
# require(httr, quietly = T)
# sp3.getlink.cddis <- function(timestamp) {
#   #internal function
#   to.GPS.timestamp <- function(year, month, day, hour, min, sec) {
#     return(as.numeric(strptime(
#       paste(year, "-", month, "-", day, " ", hour, ":", min, ":", sec, sep = ""),
#       format = "%Y-%m-%d %H:%M:%OS",
#       tz = "+0"
#     )) - 315964782 + 12)
#   }
#   
#   annee_debut = timestamp %/% (365 * 24 * 60 * 60) + 1980
#   sec_debut = timestamp
#   week_debut = sec_debut %/% (60 * 60 * 24 * 7)
#   day_of_year = (timestamp - as.numeric(to.GPS.timestamp(annee_debut, 1, 1, 0, 0, 0))) %/%
#     86400 + 2
#   temp.link = paste0(week_debut,
#                      "/", sp3_suffix,
#                      annee_debut,
#                      sprintf("%03d", day_of_year),
#                      "*")
#   
#   temp.link2 = paste0(sp3_ftp, temp.link, ".SP3.gz*?list")
#   print(paste("Trying to ping [", temp.link2, "]"))
#   result <-
#     RETRY(
#       "GET",
#       temp.link2,
#       authenticate("phuongbac", "Phuongbac2", type = "basic"),
#       times = 5,
#       quiet = F,
#       pause_base = 5
#     )
#   
#   result = (strsplit(content(result), "\r*\n")[[1]] %>% strsplit("    "))[[1]][1]
#   result = paste0(sp3_ftp, week_debut, "/", result)
#   print(paste("Got: [", result, "]"))
#   
#   return(result)
# }
# 
# sp3.download.cddis <- function(link, out.path) {
#   print("Start to download sp3 from cddis")
#   file.name = tail(strsplit(link, "/")[[1]], 1)
#   r <- RETRY(
#     "GET",
#     link,
#     write_disk(paste0(out.path, "/", file.name)),
#     authenticate("phuongbac", "Phuongbac2", type = "basic"),
#     times = 5,
#     quiet = F,
#     pause_base = 5
#   )
#   print(http_status(r)$message)
#   if (!http_error(r)) {
#     R.utils::gunzip(paste0(out.path, "/", file.name),
#                     overwrite = T,
#                     remove = T)
#     return(file.name)
#   }
#   return(NA)
# }



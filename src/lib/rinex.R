require(stringr)
require(dplyr)
require(tictoc)
#rinex = readLines(file)
rinex.getObsHeader <- function(rinexFile) {
  lines <- str_trim(readLines(rinexFile, 100))
  lines <- lines[1:which (lines == "END OF HEADER")]
  nm <- trimws(substr(lines, 61, 80))
  vv <- trimws(substr(lines, 1, 60))
  vv <- stats::setNames(vv, nm)
  return(vv)
}

rinex.getObsData <- function(rinexFile, header) {
  lines = str_trim(readLines(rinexFile))
  lines = lines[length(header) + 1:length(lines)]
  return(lines)
}

rinex.getVersion <- function(header) {
  return(as.numeric(strsplit(header["VERSION / TYPE"], " ")[[1]][1]))
}

rinex.getObsInfo <- function(header) {
  obs.info = strsplit(gsub("\\s+", " ", header[names(header) == "SYS / # / OBS TYPES"]), " ")
  constel = lapply(obs.info, "[[", 1)
  obs = lapply(
    obs.info,
    FUN = function(x) {
      return(x[3:(as.numeric(x[2]) + 2)])
    }
  )
  names(obs) = constel
  return(obs)
}

rinex.parseObs <- function(obs.data, obs.type) {
  max.len = max(lengths(obs.type)) + 1
  res = obs.data[lapply(X = obs.data, substr , 1, 1) %in% names(obs.type)] %>%
    stringr::str_squish() %>%
    strsplit(" ")
  res = lapply(res, "length<-", max.len) %>%
    unlist() %>%
    matrix(7) %>%
    t() %>% data.frame()
  time = obs.data[which(lapply(obs.data, substr, 1, 1) == ">")]
  times = substr(time, 34, 35)
  time = substr(time , 3, 21) %>%
    as.POSIXct(format = "%Y %m %d %H %M %S") - gps_time_offset
  timestamp = as.numeric(rep(time, times))
  
  res$timestamp = timestamp
  res$X1 = paste0("P", res$X1)
  return(list(data = res, obs.type = obs.type))
}

rinex.parseObs2 <- function(rinexFile) {
  # read rinex file
  obs.header = rinex.getObsHeader(rinexFile)
  obs.type = rinex.getObsInfo(obs.header)
  obs.data = rinex.getObsData(rinexFile, obs.header)
  
  # remove space in ID
  take = substr(obs.data, 2, 2) == " "
  take = take[-which(is.na(take))]
  obs.data[take] = sub("(.{1}).", "\\10", obs.data[take])
  
  # ...
  ID = lapply(obs.data, substr, 1, 3)
  max.len = max(lengths(obs.type)) + 1
  res = obs.data[lapply(X = obs.data, substr , 1, 1) %in% names(obs.type)] %>%
    stringr::str_squish() %>%
    strsplit(" ")
  res = lapply(res, "length<-", max.len) %>%
    unlist() %>%
    matrix(max.len) %>%
    t() %>% data.frame()
  time = obs.data[which(lapply(obs.data, substr, 1, 1) == ">")]
  times = substr(time, 34, 35)
  time = substr(time , 3, 21) %>%
    as.POSIXct(format = "%Y %m %d %H %M %S", tz = "+0") - 315964782
  timestamp = as.numeric(rep(time, times))
  
  res$TimeStamp = timestamp
  res$X1 = paste0("P", res$X1)
  names(res)[1]  = "ID"
  return(list(data = res, obs.type = obs.type))
}

RINEX.get.header <- function(dat) {
  for (i in 1:length(dat)) {
    line = dat[i]
    if (length(grep("END OF HEADER", line)) == 1) {
      return(dat[1:i])
    }
  }
}

RINEX.get.obs.data <- function(lines, obs.type) {
  lines_date <-
    strsplit(gsub("\\s+", 
                  " ",
                  gsub(
                    "^\\s+|\\s+$",
                    "",
                    lines[substr(lines, 1, 1) == '>']
                    )
                  ), ' ')
  annee <- sapply(lines_date, '[[', 2)
  mois <- sapply(lines_date, '[[', 3)
  jour <- sapply(lines_date, '[[', 4)
  heure <- sapply(lines_date, '[[', 5)
  minute <- sapply(lines_date, '[[', 6)
  seconde <- sapply(lines_date, '[[', 7)
  options(digits.sec = 2)
  ts <-
    as.double(strptime(
      paste(
        annee,
        '-',
        mois,
        '-',
        jour,
        ' ',
        heure,
        ':',
        minute,
        ':',
        seconde,
        sep = ''
      ),
      format = "%Y-%m-%d %H:%M:%OS",
      tz = 'UTC'
    )) - 315964800
  num.sat <- as.numeric(sapply(lines_date, '[[', 9))
  ts = rep(ts, num.sat)
  
  obs = sapply(obs.type, "[[", 1)
  line.data = lines[substr(lines, 1, 1) %in% obs]
  constel = substr(line.data,1,1)
  return(data.frame("data" = line.data, "timestamp" = ts, "constel" = constel))
}


require(dplyr)
require(tidyr)
path_input  = "input/2143-4.obs"
dat = readLines(path_input)
header = RINEX.get.header(dat)
nm <- gsub(" ", "", substr(header, 61, 80))

organization = substr(header[which(nm == "SYS/#/OBSTYPES")], 1, 60)
obs.type = strsplit(gsub("\\s+", " ", organization), " ")

obs.data = RINEX.get.obs.data(dat[length(header):length(dat)], obs.type) %>% group_by(constel) %>% nest()

RINEX = IRanges::DataFrameList()
for (i in 1:length(obs.data$constel)) {
  obs = obs.type[lapply(obs.type, "[[", 1) == obs.data$constel[i]][[1]]
  obs[2] = "PRN"
  data = strsplit(gsub("\\s+", " ", obs.data$data[[i]]$data), " ")
  lapply(data, FUN = function(x){
    x[2:length(x)] = as.numeric(x[2:length(x)])
    names(x) = obs[2:length(obs)]
    x
  })
}

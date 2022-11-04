source("src/app/gnssR-cfg.R")

if (do_extract) {
  ##############################################################################
  ######################      PREPARE FOR EXTRACT      #########################
  ##############################################################################
  print("cleaning up...")
  if (!dir.exists(extract_path))
    dir.create(extract_path)
  #do.call(file.remove, list(list.files(extract_path, full.names = TRUE)))
  rm(list = ls())
  gc()
  extract.time = Sys.time()
  ##############################################################################
  ##############################################################################
  print("--- start to extract ---")
  require(fst)
  source("src/app/gnssR-cfg.R")
  source("src/lib/rinex.R")
  source("src/lib/sp3.R")
  source("src/func/calculate_azimuth_elevation.R")
  
  list.rinex = list.files(input_path, pattern = ".obs")
  for (rnx.file in list.rinex) {
    print(paste("--- reading ephemeris from [", rnx.file, "] ---"))
    rnx = rinex.parseObs2(rinexFile = paste0(input_path, rnx.file))
    
    time = max(rnx$data$TimeStamp)
    print(paste(
      "--- Downloading corressponding sp3 file with timestamp [",
      time,
      "] ---"
    ))
    sp3.link = sp3.getLink(time)
    sp3.file = paste0(sp3_path, sp3.download(sp3.link, sp3_path))
    sp3.dat = SP3.spline(sp3.file)
    sat.data <- filter(sp3.dat,
                       TimeStamp >= (min(rnx$data$TimeStamp) - 50) &
                         TimeStamp <= (max(rnx$data$TimeStamp) + 50))
    data.merge = merge(x = sat.data,
                       y = rnx$data,
                       all = F)
    print("--- calculating elevation and azimuth for each sat ---")
    data.merge = tibble::add_column(
      data.merge,
      calculate_azimuth_elevation(data.merge$X,
                                  data.merge$Y,
                                  data.merge$Z,
                                  rec["x"],
                                  rec["y"],
                                  rec["z"])
    )
    
    rm(sat.data, sp3.dat)
    constel = substring(data.merge$ID, 2, 2)
    print("--- exporting... ---")
    for (name in names(rnx$obs.type)) {
      tmp.df = data.merge[constel == name,]
      names(tmp.df)[6:(5 + lengths(rnx$obs.type[name]))] <-
        rnx$obs.type[name][[1]]
      for (x in unique(tmp.df$ID)) {
        dplyr::filter(tmp.df, ID == x) %>%
          write.fst(paste0(
            extract_path,
            strsplit(rnx.file, ".obs")[[1]],
            "_",
            x,
            "RNX.fst"
          ))
      }
    }
  }
  print("--- finished ---")
  print(Sys.time() - extract.time)
}

if (do_analysis) {
  ##############################################################################
  ######################      PREPARE FOR ANALYSIS      ########################
  ##############################################################################
  print("cleaning up...")
  
  # do.call(file.remove, list(list.files(paste0(result_path, "tmp/"), full.names = TRUE)))
  # do.call(file.remove, list(list.files(paste0(result_path, "tmp/lsp/"), full.names = TRUE)))
  if (!dir.exists(result_path))
    dir.create(result_path)
  if (!dir.exists(extract_backup_path))
    dir.create(extract_backup_path)
  dir.create(paste0(result_path, "tmp/lsp/"), recursive=T, showWarnings=F)
  rm(list = ls())
  gc()
  analysis.time = Sys.time()
  ##############################################################################
  source("src/app/gnssR-cfg.R")
  source("src/func/determine_f_rinex.R")
  require(fst)
  require(dplyr)
  require(dplyr)
  require(data.table)
  require(pbmcapply)
  list.extract = list.files(extract_path, pattern = "RNX.fst$")
  if (length(list.extract) == 0) {
    file.copy(list.files(paste0(extract_path, "backup/"), full.names = TRUE), extract_path)
    list.extract = list.files(extract_path, pattern = "RNX.fst")
  }
  print(paste("--- start to analysis [", length(list.extract), "] ---"))

  final = pbmclapply(1:length(list.extract), FUN =function(analysis.index) {
      ephemeris.data = as.data.frame(fst(paste0(extract_path, list.extract[analysis.index])))
      if (nrow(ephemeris.data) < 3)
        return(NA)

      ephemeris.data$elevation_point = c(
        diff(ephemeris.data$elevation),
        ephemeris.data$elevation[nrow(ephemeris.data)] -
          ephemeris.data$elevation[nrow(ephemeris.data) -
                                     1]
      )
      # combine L1, L2
      signal.col = names(ephemeris.data)[which(names(ephemeris.data) %like% "^S")]
      if (length(signal.col) >= 2) {
        SXX = 0
        for (s in signal.col) {
          SXX = SXX + 10 ^ (as.numeric(ephemeris.data[, s]) / 10)
        }
        ephemeris.data$S9X = 10 * log10(SXX)
      }
      # end combine
      determine_f(ephemeris.data)
  }, mc.cores = 8)
  # ##############################################################################
  # for (analysis.index  in 1:length(list.extract)) {
  #   print(paste(analysis.index, list.extract[analysis.index]))
  #   ephemeris.data = as.data.frame(fst(paste0(extract_path, list.extract[analysis.index])))
  #   if (nrow(ephemeris.data) < 3)
  #     return(NA)
  # 
  #   ephemeris.data$elevation_point = c(
  #     diff(ephemeris.data$elevation),
  #     ephemeris.data$elevation[nrow(ephemeris.data)] -
  #       ephemeris.data$elevation[nrow(ephemeris.data) -
  #                                  1]
  #   )
  #   # combine L1, L2
  #   signal.col = names(ephemeris.data)[which(names(ephemeris.data) %like% "^S")]
  #   if (length(signal.col) >= 2) {
  #     SXX = 0
  #     for (s in signal.col) {
  #       SXX = SXX + 10 ^ (as.numeric(ephemeris.data[, s]) / 10)
  #     }
  #     ephemeris.data$S9X = 10 * log10(SXX)
  #   }
  #   # end combine
  #   determine_f(ephemeris.data)
  # }
  # ##############################################################################
  
  
  
  final = final[!is.na(final)] %>% bind_rows()
  
  write.table(final,
              file = paste0(result_path, "RESULT.txt"),
              row.names = F)
  print("--- done --- ")
  print(Sys.time() - analysis.time)
}

if (F) {
  x = filter(final, frequency != "NaN")
  x$height = as.numeric(x$frequency) * as.numeric(x$lambda) / 2
  x$time = as_datetime(as.numeric(x$timestamp), origin = '1980-1-6')
  ggplot(x) + geom_line(aes(x = time, y = height, col = signal)) + geom_point(aes(x = time, y = height, col = signal))
}
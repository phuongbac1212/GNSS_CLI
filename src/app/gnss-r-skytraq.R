source("src/app/gnssR-cfg.R")                                                                          

if (do_extract) {
  ##############################################################################
  ######################      PREPARE FOR EXTRACT      #########################
  ##############################################################################
  print("cleaning up...")
  if (!dir.exists(extract_path))
    dir.create(extract_path)
  #do.call(file.remove, list(list.files(extract_path, full.names = TRUE)))
  if (!dir.exists(input_backup_path)) 
    dir.create(input_backup_path)
  rm(list = ls())
  gc()
  extract.time = Sys.time()
  ##############################################################################
  
  print("--- start to extract ---")
  
  require(fst)
  fst::threads_fst(NULL, F)
  
  source("src/app/gnssR-cfg.R")
  source("src/lib/skytraq.R")
  source("src/lib/sp3.R")
  source("src/func/calculate_azimuth_elevation.R")
  
  list.extract = list.files(input_backup_path, pattern = ".dat$")
  if (length(list.extract) == 0) {
    file.copy(list.files(input_backup_path, full.names = TRUE),
              input_path)
    list.extract = list.files(input_path, pattern = ".dat")
  }
  
  
  
  list.stq = list.files(input_path, pattern = ".dat", full.names = T)
  for (i in 1:length(list.stq)) {
    stq.file = list.stq[i]
    date = strsplit(strsplit(stq.file, "/")[[1]][3], ".dat")[[1]][1]
    print(paste("reading ephemeris from", stq.file))
    
    stq = skytraq.parseAll(path = stq.file, type = 0xe5) %>% data.frame()
    
    time = max(stq$Timestamp)
    print(paste("Downloading corressponding sp3 file with week [", time, "]"))
    sp3.link = sp3.getLink(time)
    sp3.file = paste0(sp3_path, sp3.download(sp3.link, sp3_path))
    sp3.dat = SP3.spline(sp3.file)
    
    print("--- merging with satellite orbit ---")
    data.merge = inner_join(stq,
                            sp3.dat,
                            by = c("Timestamp" = "TimeStamp", "Observation.Type" =
                                     "ID"))
    
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
    
    
    print("--- exporting... ---")
    for (name in unique(data.merge$`Observation.Type`)) {
      tmp.df = data.merge[data.merge$`Observation.Type` == name, ]
      for (x in unique(tmp.df$`Signal.Type`)) {
        dplyr::filter(tmp.df, `Signal.Type` == x) %>%
          write.fst(paste0(extract_path, date, "_", name, "_", x, ".fst"))
      }
    }
    
    fs::file_move(stq.file, input_backup_path)
    rm(data.merge, sp3.dat, stq, tmp.df)
    closeAllConnections()
    gc()
  }
  print("--- finished ---")
  print(Sys.time() - extract.time)
}



if (do_analysis) {
  ##############################################################################
  ######################      PREPARE FOR ANALYSIS      ########################
  ##############################################################################
  print("cleaning up...")
  if (!dir.exists(result_path))
    dir.create(result_path)
  if (!dir.exists(extract_backup_path))
    dir.create(extract_backup_path)
  do.call(file.remove, list(list.files(
    paste(result_path, "tmp/"), full.names = TRUE
  )))
  dir.create(paste0(result_path, "tmp"), showWarnings = F)
  rm(list = ls())
  gc()
  analysis.time = Sys.time()
  ##############################################################################
  source("src/app/gnssR-cfg.R")
  source("src/func/determine_f.R")
  source("src/func/pos2geo.R")
  source("src/func/location_filter.R")
  source("src/func/export.R")
  source("src/func/fresnel_filter.R")
  require(foreach)
  require(fst)
  require(pbmcapply)
  require(pbapply)
  require(dplyr)
  require(data.table)
  
  list.extract = list.files(extract_path, pattern = ".fst", full.names = T)
  print(paste("--- start to analysis [", length(list.extract), "] ---"))
  
  # for (file in list.extract) {
  #   ephemeris.data = as.data.frame(fst(file))
  #   if (nrow(ephemeris.data) == 0) return(NA)
  #   elevation_point = diff(ephemeris.data$elevation)
  #   ephemeris.data = ephemeris.data[1:(nrow(ephemeris.data)-1), ]
  #   ephemeris.data$elevation_point = elevation_point
  #   s = NA#signal
  #   print(determine_f(ephemeris.data) %>% bind_rows())
  # }
  x = pbmclapply(
    list.extract,
    FUN = function(file) {
      ephemeris.data = as.data.frame(fst(file))
      if (nrow(ephemeris.data) < 3)
        return(data.frame())
      elevation_point = diff(ephemeris.data$elevation)
      ephemeris.data = ephemeris.data[1:(nrow(ephemeris.data) - 1),]
      ephemeris.data$elevation_point = elevation_point
      s = NA#signal
      return(determine_f(ephemeris.data) %>% bind_rows())
    },
    mc.cores = n.cluster
  ) %>% rbindlist()
  
  write.csv(x,
            file = paste0(result_path, "Result.csv"), row.names = F)
  print(paste("--- finished [", length(list.files(
    paste0(result_path, "tmp/")
  )), "] file created ---"))
  print(Sys.time() - analysis.time)
  
  system(paste0("gedit ", result_path, "Result.csv"))
}

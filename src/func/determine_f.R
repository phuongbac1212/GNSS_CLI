source("src/func/calcul_lsp.R")
source("src/func/cal_dlat_dlon.R")

determine_f_min_max <-
  function(elevation_cur,
           elevation_point_cur,
           longueur_onde) {
    f_min <- 99999999999999
    f_max <- 0
    f_tot <- 0
    h_point_min <- -h_point_max
    
    f1 <-
      2 / longueur_onde * (h_min + h_point_max / elevation_point_cur * tan(elevation_cur))
    f2 <-
      2 / longueur_onde * (h_max + h_point_max / elevation_point_cur * tan(elevation_cur))
    f3 <-
      2 / longueur_onde * (h_min - h_point_max / elevation_point_cur * tan(elevation_cur))
    f4 <-
      2 / longueur_onde * (h_max - h_point_max / elevation_point_cur * tan(elevation_cur))
    
    f_max <- max(abs(c(f1, f2, f3, f4)))
    f_min <- min(abs(c(f1, f2, f3, f4)))
    
    if (f1 < 0 & f2 < 0 & f3 < 0 &  f4 < 0) {
      signe <- -1
    }
    else{
      signe <- 1
    }
    
    
    return(c(f_min, f_max, signe))
  }


calcul_f_a_phi <-
  function(ma_date,
           snr,
           elevation,
           elevation_point,
           longueur_onde,
           f_min,
           f_max,
           longueur_minimum_serie) {
    # Detrendage du snr  --> on enleve la contribution du signal direct.
    if (length(snr) > longueur_minimum_serie) {
      snr_detrended <- detrend_snr(ma_date, snr)
    }
    else{
      snr_detrended <- snr
    }
    f <- NaN
    amplitude <- NaN
    phase <- NaN
    f <-
      calcul_temps_frequence(sin(elevation), snr_detrended, f_min, f_max, 'lsp') #methode possible : 'lsp' ou 'RIAA'
    
    if (h_fixed) {
      f_calcul <- 2 * h0 / longueur_onde
    }
    else{
      f_calcul <- f
    }
    # Determination de l'amplitude et du dephasage par les moindres carres
    K1 <- cos(2 * pi * f_calcul * sin(elevation))
    K2 <- -sin(2 * pi * f_calcul * sin(elevation))
    K3 <- snr_detrended
    P <- diag(length(K1))
    matA <- mat.or.vec(length(K1), 2)
    matA[, 1] <- K1
    matA[, 2] <- K2
    mamat <- NA
    require(MASS)
    try(mamat <- ginv(t(matA) %*% P %*% matA), silent = TRUE)
    if (length(mamat) > 1) {
      result <- mamat %*% (t(matA) %*% P %*% K3)
      Acos_phi <- result[1]
      Asin_phi <- result[2]
      
      v <- matA %*% c(result[1], result[2]) - K3
      emq_0 <- sqrt((t(v) %*% P %*% v) / (length(K1) - 2))
      amplitude <- sqrt(Acos_phi ^ 2 + Asin_phi ^ 2)
      
      acos_phi <- acos(Acos_phi / amplitude)
      asin_phi <- asin(Asin_phi / amplitude)
      if (asin_phi < 0) {
        phase <- -acos_phi
      }
      else if (asin_phi > 0) {
        phase <- acos_phi
      }
      phase <- phase %% (2 * pi)
      
    }
    else{
      phase <- NaN
      amplitude <- NaN
    }
    return(c(amplitude, phase, f))
  }

determine_f <- function(ephemeris.data) {
  source('src/func/detrend_snr.R')
  source('src/func/decoupe.R')
  
  require(MASS, quietly = T)
  require(lubridate, quietly = T)
  require(tidyr)
  
  pnr = unique(ephemeris.data$Observation.Type)
  # on apply filter on elevation
  ephemeris.data = subset(ephemeris.data,
                          elevation < elevation_max &
                            elevation > elevation_min)
  
  # on apply filer on azimuth
  temp.index = 0
  for (i in seq(1, length(azimuth_to_take), 2)) {
    temp.index <-
      c(
        temp.index,
        which(
          ephemeris.data$azimuth >= azimuth_to_take[i] &
            ephemeris.data$azimuth <= azimuth_to_take[i + 1]
        )
      )
  }
  ephemeris.data <- ephemeris.data[unique(temp.index), ]
  rm(temp.index)
  
  # determine f
  ephemeris.data$snr = 10 ^ (as.numeric(ephemeris.data$CNR) / 20) # convert dbHz to Volt
  
  
  a_prendre <-
    decoupe(
      ephemeris.data$Timestamp,
      ephemeris.data$elevation,
      TRUE,
      max_hole,
      minimum_length_series
    )
  #if(is.na(a_prendre)) return(NA)
  if (length(a_prendre) != 0) {
    res = lapply(
      a_prendre,
      FUN = function(take) {
        if(length(take)==1) 
          if(is.na(take)) 
            return(NA)
        determine.f.data = ephemeris.data[take,]
        longueur_onde = as.numeric(unique(determine.f.data$Wavelength))
        if (length(longueur_onde) > 1) {
          return(NA)
        }
        if (method == 'Larson') {
          f_min <- 2 / longueur_onde * h_min
          f_max <- 2 / longueur_onde * h_max
          s <<- unique(determine.f.data$`Signal Type`)
          result <-
            calcul_f_a_phi(
              determine.f.data$Timestamp,
              determine.f.data$snr,
              determine.f.data$elevation,
              determine.f.data$elevation_point,
              longueur_onde,
              f_min,
              f_max,
              minimum_length_series
            )
          result2 <-
            list(
              "timestamp" = as.numeric(mean(as.numeric(
                determine.f.data$Timestamp
              ))),
              "PRN" = pnr,
              "elevation" = mean(determine.f.data$elevation),
              "elevation_point" = mean(determine.f.data$elevation_point),
              "azimuth" = mean(determine.f.data$azimuth),
              "amplitude" = result[1],
              "phase" = result[2],
              "wavelength" = longueur_onde,
              "frequency" = result[3],
              "height" = result[3] * longueur_onde / 2,
              "signal" = s
            )
          # result2 = append(result2, cal.fresnel.ellipse(h = result[3] *longueur_onde /2, mean(elevation_cur), longueur_onde))
          # result2 = append(result2, cal.dlat.dlon(pos_E_xyz= ,elevation = mean(elevation_cur), azimuth = mean(azimuth_cur))
          
          result2 = append(
            result2,
            cal.dlat.dlon(
              pos_E_xyz = c(mean(determine.f.data$X), mean(determine.f.data$Y), mean(determine.f.data$Z)),
              pos_R_xyz = rec,
              elevation = mean(determine.f.data$elevation),
              azimuth = mean(determine.f.data$azimuth),
              wavelength = longueur_onde
            )
          )
          names(result2)[12:26] <-
            c(
              "Lon_simu",
              "Lat_simu",
              "Height_simu",
              "Fz_a",
              "Fz_b",
              "Fz_c",
              "A_x",
              "A_y",
              "A_z",
              "B_x",
              "B_y",
              "B_z",
              "AB_x",
              "AB_y",
              "AB_z",
              "signal"
            )
          
          return(result2)
        }
      }
    )# end lapply
    return (res[!is.na(res)])
  }
}# end func

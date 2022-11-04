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
    library(MASS)
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
    return(c(
      "amplitude" = amplitude,
      "phase" = phase,
      "frequency" = f
    ))
  }

determine_f <- function(ephemeris.data) {
  source('src/func/detrend_snr.R')
  source('src/func/decoupe.R')
  source('src/func/determine_lambda.R')
  source('src/func/detrend_snr.R')
  source('src/func/calcul_lsp.R')
  source("src/func/cal_dlat_dlon.R")
  require(MASS, quietly = T)
  require(lubridate, quietly = T)
  require(tidyr)
  require(dplyr)
  require(foreach)
  
  prn = unique(ephemeris.data$ID)
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
  ephemeris.data <- ephemeris.data[unique(temp.index),]
  rm(temp.index)
  
  if (nrow(ephemeris.data) < minimum_length_series) {
    paste(paste("   +++ not enough data left after filter elevation and azimuth ", prn, "+++   "))
    return(NA)
  }
  list_f = decoupe(
    ephemeris.data$TimeStamp,
    ephemeris.data$elevation_point,
    T,
    max_hole,
    minimum_length_series
  )
  if (length(list_f) == 0)
    return(NA)
  else if (length(list_f) == 1)
    if (is.na(list_f)) return(NA)
  
  res = lapply(
    list_f,
    FUN = function(i) {
      ep.temp = ephemeris.data[i, ]
      signal.names = names(ep.temp)[which(names(ephemeris.data) %like% "^S[0-9]")]
      total.result = NA
      for (signal.name in signal.names) {
        print(signal.name)
        if (sum(is.na(ep.temp[, signal.name])) > max_hole) {
          return(NA)
        } else if (sum(is.na(ep.temp[, signal.name])) > 0) {
          ep.temp = drop_na(ep.temp)
        }
        
        snr = 10 ^ (as.numeric(ep.temp[, signal.name]) / 20)
        date =  ep.temp$TimeStamp
        elevation <- as.numeric(ep.temp$elevation)
        azimuth <- as.numeric(ep.temp$azimuth)
        elevation_point <- ep.temp$elevation_point
        # elevation_point <-
        #   c(elevation_point, elevation_point[length(elevation_point)]) # add final e.
        
        if (method == "Larson") {
          lambda = determine_lambda(prn, signal.name, prn, file.type = "rinex")
          
          f_min <- 2 / lambda * h_min
          f_max <- 2 / lambda * h_max
          result <-
            calcul_f_a_phi(
              date,
              snr,
              elevation,
              elevation_point,
              lambda,
              f_min,
              f_max,
              minimum_length_series
            )
          dlat_dlon = cal.dlat.dlon(c(mean(ephemeris.data$X), mean(ephemeris.data$Y), mean(ephemeris.data$Z)), 
                                    rec, 
                                    elevation = mean(elevation_cur),
                                    azimuth = mean(azimuth_cur),
                                    wavelength = lambda)

          total.result = append(total.result, list(
            c(
              timestamp = mean(date),
              prn = prn,
              signal = signal.name,
              lambda = lambda,
              elevation = mean(elevation),
              elevation_dot = mean(elevation_point),
              azimuth = mean(azimuth),
              result,
              dlat_dlon
            )
          ))
        }
      }
      return(total.result)
    }
  )[[1]]
  if(length(res) == 1)
    if (is.na(res)) return(res)
  return(res[!is.na(res)] %>% bind_rows())
}

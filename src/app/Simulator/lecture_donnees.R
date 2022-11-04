lecture_donnees <- function (nom_sat_liste,mon_data,ma_date,organisation){
      data_final <- list()
      nom_snr <- 0

      for(u in 1:length(nom_sat_liste)){
            data_final[[u]] <- list()
            condition <- gsub(' ','0',mon_data[,1])==nom_sat_liste[u] #sapply(mon_data,'[[',1) == nom_sat_liste[u]
            data_final[[u]][[1]] <- nom_sat_liste[u]
            data_final[[u]][[2]] <- ma_date[condition]
            incr <- 3
            for(i in 1:length(organisation)){
                  nom_snr_cur <- organisation[i]
                  if(nom_snr_cur %in% liste_codes_possibles & sum(condition)!=0){
                        mon_data <<- mon_data
                        organisation <<- organisation
                        SNR_cur <- as.numeric(mon_data[condition,][,i-1])#sapply(mon_data[condition][sapply(mon_data[condition],length)>=i-1],'[[',i-1)  
                        date_cur <- ma_date[condition]
                        SNR_cur[is.na(SNR_cur)] <- 1 
                        if(length(date_cur)==0){
                              SNR_cur <- as.numeric(data_final[[u]][[incr-1]])*0+1
                        }
                        else{
                              SNR_cur[SNR_cur==''] <- 1
                              SNR_cur <- approx(date_cur,SNR_cur,data_final[[u]][[2]],method='constant',rule=1)$y
                        }
                              data_final[[u]][[incr]] <- SNR_cur

                        nom_snr[incr-2] <- nom_snr_cur
                        incr <- incr + 1
                  }
            }
            data_final[[u]][[incr]] <- nom_snr
      }

      return(data_final)
}
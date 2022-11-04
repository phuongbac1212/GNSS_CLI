calcul_f_a_phi <- function(ma_date,snr,elevation,elevation_point,longueur_onde,f_min,f_max,longueur_minimum_serie){
      # Detrendage du snr  --> on enleve la contribution du signal direct.
      if(length(snr)>longueur_minimum_serie){
            snr_detrended <- detrend_snr(ma_date,snr)
      }
      else{snr_detrended <- snr}
      f <- NaN
      amplitude <- NaN
      phase <- NaN          
      f <- calcul_temps_frequence(sin(elevation),snr_detrended,f_min,f_max,'lsp') #methode possible : 'lsp' ou 'RIAA'
          
      if(h_fixed){
            f_calcul <- 2*H0/longueur_onde
      }
      else{
            f_calcul <- f
      }
      # Determination de l'amplitude et du dephasage par les moindres carres
      K1 <- cos(2*pi*f_calcul*sin(elevation))
      K2 <- -sin(2*pi*f_calcul*sin(elevation))
      K3 <- snr_detrended
      P <- diag(length(K1))
      matA <- mat.or.vec(length(K1),2)
      matA[,1] <- K1
      matA[,2] <- K2
      mamat <- NA
      library(MASS)
      try(mamat <- ginv(t(matA)%*%P%*%matA),silent=TRUE)   
      if(length(mamat)>1){
            result <- mamat%*%(t(matA)%*%P%*%K3)
            Acos_phi <- result[1]
            Asin_phi <- result[2]

            v <- matA%*%c(result[1],result[2])-K3
            emq_0 <- sqrt((t(v)%*%P%*%v)/(length(K1)-2))
            amplitude <- sqrt(Acos_phi^2+Asin_phi^2)

            acos_phi <- acos(Acos_phi/amplitude)
            asin_phi <- asin(Asin_phi/amplitude)
            if(asin_phi <0){phase <- -acos_phi}
            else if(asin_phi >0){phase <- acos_phi}
            phase <- phase %% (2*pi)
            
      }
      else{phase <- NaN;amplitude <- NaN}
      return(c(amplitude,phase,f)) 
}

determine_f <- function(nom_fichier,autonome){

      # Importation des fonctions necessaires
      source('detrend_snr.r')
      source('determine_h.r')
      source('determine_lambda.r')
      source('determine_f_min_max.r') 
      source('decoupe.r')
      source('calcul_lsp.r')
      library(MASS,quietly=TRUE)
      library(lubridate,quietly=TRUE)

      if(autonome){
            source('constants.r') 
            source('lecture_directeur.r') 
      }
        
      date_final <- 0
      nom_sat_final <- 0
      elevation_final <- 0
      azimuth_final <- 0
      elevation_point_final <- 0
      amplitude_final <- 0
      phase_final <- 0
      longueur_onde_final <- 0
      f_final <- 0
      incr <- 1
      
      organisation <- scan(paste(path_extraction,nom_fichier,sep=''),what=character(0),nlines=1,quiet=TRUE,skip=0)
      nb_col_organisation <- length(organisation)
      nom_sat <- substr(nom_fichier,nchar(nom_fichier)-6,nchar(nom_fichier)-4)
      if(nom_sat %in% nom_sat_liste){
            mon_data <- matrix(scan(paste(path_extraction,nom_fichier,sep=''),what=character(0),skip=1,quiet=TRUE),ncol=nb_col_organisation,byrow=TRUE)
   			    if(autonome){
        					mon_data[,8] <- as.numeric(mon_data[,8])*pi/180
        					mon_data[,9] <- as.numeric(mon_data[,9])*pi/180
        					mon_data[,11][as.numeric(mon_data[,11])==0] <- '-----' 
        					mon_data[,12][as.numeric(mon_data[,12])==0] <- '-----' 
        					mon_data[,13][as.numeric(mon_data[,13])==0] <- '-----' 
        					mon_data[,14][as.numeric(mon_data[,14])==0] <- '-----' 
		        }
            # on applique les filtres en elevations

            mon_data <- mon_data[(as.numeric(mon_data[,8])>elevation_min & as.numeric(mon_data[,8])<elevation_max),]
            mon_data <- mon_data[!is.na(mon_data[,1]),]
            # on applique les filtres en azimuth
            if(dim(mon_data)[1]!=0){
                  boolean <- seq(1,dim(mon_data)[1],1)*0  
                  seq_a_prendre <- seq(1,dim(mon_data)[1],1)    
                  azi <- as.numeric(mon_data[,9]) %% (2*pi)
                  for(u in seq(1,length(masque_azimuth),2)){
                        boolean[(azi>=masque_azimuth[u] & azi <= masque_azimuth[u+1])] <- boolean[(azi>=masque_azimuth[u] & azi <= masque_azimuth[u+1])]+1
                  }
                  boolean[boolean!=0] <- 1
                  seq_a_prendre <- seq_a_prendre[as.logical(boolean)]
                  mon_data <- mon_data[seq_a_prendre,]
                 
                  ##########  Determination de f
                  
                  print(paste('Analyse du fichier ',nom_fichier,' avec la methode ',methode_a_utiliser,sep=''),quote=FALSE)
      
                  for(num_code in 1:length(liste_code_a_utiliser)){
                        code_a_utiliser <<- liste_code_a_utiliser[num_code]
                        if(length(mon_data)==length(organisation)){snr <- as.numeric(mon_data[organisation==code_a_utiliser])}
                        else{snr <- as.numeric(mon_data[,organisation==code_a_utiliser])}
                        if(length(snr)>longueur_minimum_serie){
                              seq_a_prendre <- seq(1,length(snr),1)
                              seq_a_prendre <- seq_a_prendre[!is.na(snr)]   
                              if(length(seq_a_prendre)!=0){
                                    snr <- 10^(snr[seq_a_prendre]/20)   # On convertit les dBHZ en Volt
                                    ma_date <- as.double(as.POSIXct(paste(mon_data[seq_a_prendre,1],'-',mon_data[seq_a_prendre,2],'-',mon_data[seq_a_prendre,3],' ',mon_data[seq_a_prendre,4],':',mon_data[seq_a_prendre,5],':',mon_data[seq_a_prendre,6],sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'),tz='UTC')
                                    if(autonome){ma_date <- ma_date - 17}  #decalage en secondes entre temps GPS et temps UTC
                                    elevation <- as.numeric(mon_data[seq_a_prendre,8])
                                    azimuth <- as.numeric(mon_data[seq_a_prendre,9]) %% (2*pi)
                                    longueur_onde <- determine_lambda(nom_sat,code_a_utiliser)
                                    elevation <- smooth.spline(ma_date,elevation)$y
                                    elevation_point <- diff(elevation)
                                    ma_date <- ma_date[1:(length(ma_date)-1)]
                                    elevation <- elevation[1:(length(elevation)-1)]
                                    azimuth <- azimuth[1:(length(azimuth)-1)]
                                    snr <- snr[1:(length(snr)-1)]
                                    a_prendre <- list()
                                    a_prendre <- decoupe(ma_date,elevation,TRUE,trou_max,longueur_minimum_serie)
                                    if(length(a_prendre)!=0){              
                                          for(u in 1:length(a_prendre)){
                                                ma_date_cur <- ma_date[a_prendre[[u]]]
                                                elevation_cur <- elevation[a_prendre[[u]]]
                                                elevation_point_cur <- elevation_point[a_prendre[[u]]]
                                                snr_cur <- snr[a_prendre[[u]]]
                                                azimuth_cur <- azimuth[a_prendre[[u]]]
                                                if(methode_a_utiliser=='LSM'){
                                                      # On decoupe encore plus si on utilise la methode LSM
                                                      ma_date_centrale <- ma_date_cur[1]
                                                      maxi <- max(ma_date_cur)
                                                      if(is.na(maxi)){maxi <- 0}
                                                      while(ma_date_centrale < maxi & !is.na(ma_date_centrale)){ # on parcourt toute la serie temporelle avec le pas definit par 1/Fh (specifie dans le fichier directeur).                                         
                                                            seq_complete <- seq(1,length(ma_date_cur),1)
                                                            index_central <- abs(ma_date_cur-ma_date_centrale)==min(abs(ma_date_cur-ma_date_centrale))
                                                            indice_cur <- seq_complete[index_central][1]
                                                            resultat <- determine_f_min_max(elevation_cur[indice_cur],elevation_point_cur[indice_cur],longueur_onde)
                                                            f_min <- resultat[1]
                                                            f_max <- resultat[2]
                                                            signe_f <- resultat[3]
                                                            T_max <- 1/f_min
                                                            T_min <- 1/f_max
                                                            seq_complete <- seq(1,length(snr_cur),1)
                                                            T <- n_period*T_max #en sin(elevation)
                                                            debut <- round(seq_complete[(abs(sin(elevation_cur)-(sin(elevation_cur[indice_cur])-T/2)))==(min(abs(sin(elevation_cur)-(sin(elevation_cur[indice_cur])-T/2))))],0)
                                                            fin <- round(seq_complete[(abs(sin(elevation_cur)-(sin(elevation_cur[indice_cur])+T/2)))==(min(abs(sin(elevation_cur)-(sin(elevation_cur[indice_cur])+T/2))))],0)
                                                            #index_min <- abs(ma_date_cur-(ma_date_centrale-((1/Fh)/2)))==min(abs(ma_date_cur-(ma_date_centrale-((1/Fh)/2))))
      #                                                      index_max <- abs(ma_date_cur-(ma_date_centrale+((1/Fh)/2)))==min(abs(ma_date_cur-(ma_date_centrale+((1/Fh)/2)))) 
      #                                                      debut <- seq_complete[index_min][1]
      #                                                      fin <- seq_complete[index_max][1]
                                                            tot <- c(debut,fin)
                                                            fin <- max(tot)
                                                            debut <- min(tot)
                   
                                                            condition <- (fin-debut)<=(1/Fh)
                                                            
                                                            condition_2 <- abs((indice_cur-debut)-(fin-indice_cur))<10

                                                            condition <- TRUE
                                                             condition_2 <- TRUE

                                                   #         if((fin-debut)>=min(T,length(ma_date)) & condition_2 & length(snr_cur[debut:fin])>longueur_minimum_serie & condition){
                                                            if(condition & condition_2){
                                                                  result <- calcul_f_a_phi(ma_date_cur[debut:fin],snr_cur[debut:fin],elevation_cur[debut:fin],elevation_point[debut:fin],longueur_onde,f_min,f_max,longueur_minimum_serie) 
                                                                  date_final[incr] <- mean(ma_date_cur[debut:fin])
                                                                  nom_sat_final[incr] <- nom_sat
                                                                  elevation_final[incr] <- mean(elevation_cur[debut:fin])
                                                                  azimuth_final[incr] <- mean(azimuth_cur[debut:fin])
                                                                  elevation_point_final[incr] <- mean(elevation_point_cur[debut:fin])
                                                                  amplitude_final[incr] <- result[1]
                                                                  phase_final[incr] <- result[2]
                                                                  longueur_onde_final[incr] <- longueur_onde
                                                                  f_final[incr] <- result[3]
                                                                  incr <- incr + 1           
                                                            }
                                                            ma_date_centrale <- ma_date_centrale + 1/Fh# ma_date[fin] #recouvrement de 50 % entre chaque fenetre mobile
                                                      }
                                                }
                                                if(methode_a_utiliser=='Larson'){
                                                      f_min <- 2/longueur_onde*h_min
                                                      f_max <- 2/longueur_onde*h_max                           
                                                      result <- calcul_f_a_phi(ma_date_cur,snr_cur,elevation_cur,elevation_point_cur,longueur_onde,f_min,f_max,longueur_minimum_serie)
                                                      date_final[incr] <- mean(ma_date_cur)
                                                      nom_sat_final[incr] <- nom_sat
                                                      elevation_final[incr] <- mean(elevation_cur)
                                                      azimuth_final[incr] <- mean(azimuth_cur)
                                                      elevation_point_final[incr] <- mean(elevation_point_cur)
                                                      amplitude_final[incr] <- result[1]
                                                      phase_final[incr] <- result[2]
                                                      longueur_onde_final[incr] <- longueur_onde
                                                      f_final[incr] <- result[3]
                                                      incr <- incr + 1  
                                                }
                                          }
                                    }
                              }
                        }
                  }
            }
      }

      date_final <- as.POSIXct(date_final,origin='1970-01-01',tz='UTC')
      col1 <- c('Year',format(as.numeric(year(date_final)),digits=1))
      col2 <- c('Month',month(date_final))
      col3 <- c('Day',day(date_final))
      col4 <- c('Hour',hour(date_final))
      col5 <- c('Minute',minute(date_final))
      col6 <- c('Seconde',format(second(date_final),digits=1))
      col7 <- c('PRN',nom_sat_final)
      col8 <- c('Elevation(rad)',format(elevation_final,digits=8))
      col9 <- c('Azimuth(rad)',format(azimuth_final,digits=8))
      col10 <- c('Elevation_point(rad/s)',format(elevation_point_final,digits=8))
      col11 <- c('Amplitude',format(amplitude_final,digits=4))
      col12 <- c('Phase',format(phase_final,digits=4))
      col13 <- c('Longueur_onde(m)',format(longueur_onde_final,digits=4))
      col14 <- c('f',format(f_final,digits=8))
      if(length(col1)>=2 & date_final != 0){
            fichier_sortie <- paste(path_output,paste('sortie_f_',nom_fichier,sep=''),sep='')    
            write.matrix(matrix(c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14),ncol=14), file = fichier_sortie,sep='\t')
      }
} # end of the function
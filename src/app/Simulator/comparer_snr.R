rm(list=ls())  # vide le cache
graphics.off()  # ferme toutes les figures deja ouvertes
options(digits.sec=2)




#
# Determination du Work Directory et des differents chemins d'acces
#

wd <- getwd()
setwd('..')
path_main <<- paste(getwd(),'/',sep='')
path_input <<- paste(getwd(),'/0.input/',sep='')
path_extraction <<- paste(getwd(),'/1.extraction/',sep='')
path_output <<- paste(getwd(),'/2.output/',sep='')
setwd(paste(getwd(),'/programs',sep='')) # We put back the Work Directory in the Programs folder


nom_sat_liste <<- matrix(scan(paste(path_main,'liste_satellites.txt',sep=''),what=character(0),quiet=TRUE),ncol=1)[,1]
elevation_min <<- 0
elevation_max <<- 90
masque_azimuth <<- c(0,0)
liste_code_a_utiliser <<- 'S1C'
trou_max <<- 10
longueur_minimum_serie <<- 300
h_min <<- 50
h_max <<- 300
longueur_onde <<- 0.1905
nom_x <- 'date' # date ou sin
volt_par_volt <- FALSE


# Importation des fonctions necessaires

source('decoupe.r')
source('calcul_lsp.r')


date_log0 <- list()
nom_sat_log0 <- 0
elevation_log0 <- list()
azimuth_log0 <- list()
SNR_log0 <- list()
incr <- 1

liste_extraction_complet_log0 <- list.files(path = path_extraction, pattern = ("log0"))
liste_extraction_complet_log1 <- list.files(path = path_extraction, pattern = ("log1"))
liste_extraction_complet_base <- list.files(path = path_extraction, pattern = ("esperce"))

for(fichier in liste_extraction_complet_log0){
      print(paste('Lecture en cours du fichier : ',fichier,sep=''),quote=FALSE)
      organisation <- scan(paste(path_extraction,fichier,sep=''),what=character(0),nlines=1,quiet=TRUE)
      nb_col_organisation <- length(organisation)
      nom_sat <- substr(fichier,nchar(fichier)-6,nchar(fichier)-4)

      if(nom_sat %in% nom_sat_liste){
            mon_data <- matrix(scan(paste(path_extraction,fichier,sep=''),what=character(0),skip=1,quiet=TRUE),ncol=nb_col_organisation,byrow=TRUE)

            # on applique les filtres en elevations
            mon_data <- mon_data[mon_data[,8]>elevation_min & mon_data[,8]<elevation_max,]

            # on applique les filtres en azimuth
            if(dim(mon_data)[1]!=0){
                  seq_a_prendre <- seq(1,dim(mon_data)[1],1)

                  for(u in seq(1,length(masque_azimuth),2)){
                        seq_a_prendre <- seq_a_prendre[!(mon_data[,9]>=masque_azimuth[u] & mon_data[,9] <= masque_azimuth[u+1])]
                  }
                  mon_data <- mon_data[seq_a_prendre,]

                  # On convertit les dBHZ en Volt
                  if(volt_par_volt){
                        mon_data[,11:nb_col_organisation] <- 10^(as.numeric(mon_data[,11:nb_col_organisation])/20)
                  }
                  
                  
                  for(num_code in 1:length(liste_code_a_utiliser)){
                        code_a_utiliser <<- liste_code_a_utiliser[num_code]
                        snr <- mon_data[,organisation==code_a_utiliser]
                        if(length(snr)!=0){
                              seq_a_prendre <- seq(1,length(snr),1)
                              seq_a_prendre <- seq_a_prendre[!is.na(snr)]
                              snr <- snr[seq_a_prendre]
                              ma_date <- as.double(strptime(paste(mon_data[seq_a_prendre,1],'-',mon_data[seq_a_prendre,2],'-',mon_data[seq_a_prendre,3],' ',mon_data[seq_a_prendre,4],':',mon_data[seq_a_prendre,5],':',mon_data[seq_a_prendre,6],sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'),tz='UTC')
                              elevation <- as.numeric(mon_data[seq_a_prendre,8])
                              azimuth <- as.numeric(mon_data[seq_a_prendre,9])
                              ma_date <- ma_date[1:(length(ma_date)-1)]
                              elevation <- elevation[1:(length(elevation)-1)]
                              azimuth <- azimuth[1:(length(azimuth)-1)]
                              snr <- snr[1:(length(snr)-1)]
                              a_prendre <- list()
                              a_prendre <- decoupe(ma_date,0,FALSE)

                              if(length(a_prendre)!=0){
                                    for(u in 1:length(a_prendre)){
                                          ma_date_cur <- ma_date[a_prendre[[u]]]
                                          elevation_cur <- elevation[a_prendre[[u]]]
                                          snr_cur <- snr[a_prendre[[u]]]
                                          azimuth_cur <- azimuth[a_prendre[[u]]]

                                          date_log0[[incr]] <- ma_date_cur
                                          nom_sat_log0[incr] <- nom_sat
                                          elevation_log0[[incr]] <- elevation_cur
                                          azimuth_log0[[incr]] <- azimuth_cur
                                          SNR_log0[[incr]] <- snr_cur
                                          incr <- incr + 1
                                    }
                              }
                        }
                  }
            }
      }
}


date_log1 <- list()
nom_sat_log1 <- 0
elevation_log1 <- list()
azimuth_log1 <- list()
SNR_log1 <- list()
incr <- 1

for(fichier in liste_extraction_complet_log1){
      print(paste('Lecture en cours du fichier : ',fichier,sep=''),quote=FALSE)
      organisation <- scan(paste(path_extraction,fichier,sep=''),what=character(0),nlines=1,quiet=TRUE)
      nb_col_organisation <- length(organisation)
      nom_sat <- substr(fichier,nchar(fichier)-6,nchar(fichier)-4)

      if(nom_sat %in% nom_sat_liste){
            mon_data <- matrix(scan(paste(path_extraction,fichier,sep=''),what=character(0),skip=1,quiet=TRUE),ncol=nb_col_organisation,byrow=TRUE)

            # on applique les filtres en elevations
            mon_data <- mon_data[mon_data[,8]>elevation_min & mon_data[,8]<elevation_max,]

            # on applique les filtres en azimuth
            if(dim(mon_data)[1]!=0){
                  seq_a_prendre <- seq(1,dim(mon_data)[1],1)

                  for(u in seq(1,length(masque_azimuth),2)){
                        seq_a_prendre <- seq_a_prendre[!(mon_data[,9]>=masque_azimuth[u] & mon_data[,9] <= masque_azimuth[u+1])]
                  }
                  mon_data <- mon_data[seq_a_prendre,]

                  # On convertit les dBHZ en Volt
                  if(volt_par_volt){
                    mon_data[,11:nb_col_organisation] <- 10^(as.numeric(mon_data[,11:nb_col_organisation])/20)
                  }
                  
                  for(num_code in 1:length(liste_code_a_utiliser)){
                        code_a_utiliser <<- liste_code_a_utiliser[num_code]
                        snr <- mon_data[,organisation==code_a_utiliser]
                        if(length(snr)!=0){
                              seq_a_prendre <- seq(1,length(snr),1)
                              seq_a_prendre <- seq_a_prendre[!is.na(snr)]
                              snr <- snr[seq_a_prendre]
                              ma_date <- as.double(strptime(paste(mon_data[seq_a_prendre,1],'-',mon_data[seq_a_prendre,2],'-',mon_data[seq_a_prendre,3],' ',mon_data[seq_a_prendre,4],':',mon_data[seq_a_prendre,5],':',mon_data[seq_a_prendre,6],sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'),tz='UTC')
                              elevation <- as.numeric(mon_data[seq_a_prendre,8])
                              azimuth <- as.numeric(mon_data[seq_a_prendre,9])
                              ma_date <- ma_date[1:(length(ma_date)-1)]
                              elevation <- elevation[1:(length(elevation)-1)]
                              azimuth <- azimuth[1:(length(azimuth)-1)]
                              snr <- snr[1:(length(snr)-1)]
                              a_prendre <- list()
                              a_prendre <- decoupe(ma_date,0,FALSE)

                              if(length(a_prendre)!=0){
                                    for(u in 1:length(a_prendre)){
                                          ma_date_cur <- ma_date[a_prendre[[u]]]
                                          elevation_cur <- elevation[a_prendre[[u]]]
                                          snr_cur <- snr[a_prendre[[u]]]
                                          azimuth_cur <- azimuth[a_prendre[[u]]]

                                          date_log1[[incr]] <- ma_date_cur
                                          nom_sat_log1[incr] <- nom_sat
                                          elevation_log1[[incr]] <- elevation_cur
                                          azimuth_log1[[incr]] <- azimuth_cur
                                          SNR_log1[[incr]] <- snr_cur
                                          incr <- incr + 1
                                    }
                              }
                        }
                  }
            }
      }
}

date_base <- list()
nom_sat_base <- 0
elevation_base <- list()
azimuth_base <- list()
SNR_base <- list()
incr <- 1

for(fichier in liste_extraction_complet_base){
      print(paste('Lecture en cours du fichier : ',fichier,sep=''),quote=FALSE)
      organisation <- scan(paste(path_extraction,fichier,sep=''),what=character(0),nlines=1,quiet=TRUE)
      nb_col_organisation <- length(organisation)
      nom_sat <- substr(fichier,nchar(fichier)-6,nchar(fichier)-4)

      if(nom_sat %in% nom_sat_liste){
            mon_data <- matrix(scan(paste(path_extraction,fichier,sep=''),what=character(0),skip=1,quiet=TRUE),ncol=nb_col_organisation,byrow=TRUE)

            # on applique les filtres en elevations
            mon_data <- mon_data[mon_data[,8]>elevation_min & mon_data[,8]<elevation_max,]

            # on applique les filtres en azimuth
            if(dim(mon_data)[1]!=0){
                  seq_a_prendre <- seq(1,dim(mon_data)[1],1)

                  for(u in seq(1,length(masque_azimuth),2)){
                        seq_a_prendre <- seq_a_prendre[!(mon_data[,9]>=masque_azimuth[u] & mon_data[,9] <= masque_azimuth[u+1])]
                  }
                  mon_data <- mon_data[seq_a_prendre,]
                  # On convertit les dBHZ en Volt
                  if(volt_par_volt){
                        mon_data[,11:nb_col_organisation] <- 10^(as.numeric(mon_data[,11:nb_col_organisation])/20)
                  }
                  
                  for(num_code in 1:length(liste_code_a_utiliser)){
                        code_a_utiliser <<- liste_code_a_utiliser[num_code]
                        snr <- mon_data[,organisation==code_a_utiliser]
                        if(length(snr)!=0){
                              seq_a_prendre <- seq(1,length(snr),1)
                              seq_a_prendre <- seq_a_prendre[!is.na(snr)]
                              snr <- snr[seq_a_prendre]
                              ma_date <- as.double(strptime(paste(mon_data[seq_a_prendre,1],'-',mon_data[seq_a_prendre,2],'-',mon_data[seq_a_prendre,3],' ',mon_data[seq_a_prendre,4],':',mon_data[seq_a_prendre,5],':',mon_data[seq_a_prendre,6],sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'),tz='UTC')
                              elevation <- as.numeric(mon_data[seq_a_prendre,8])
                              azimuth <- as.numeric(mon_data[seq_a_prendre,9])
                              ma_date <- ma_date[1:(length(ma_date)-1)]
                              elevation <- elevation[1:(length(elevation)-1)]
                              azimuth <- azimuth[1:(length(azimuth)-1)]
                              snr <- snr[1:(length(snr)-1)]
                              a_prendre <- list()
                              a_prendre <- decoupe(ma_date,0,FALSE)

                              if(length(a_prendre)!=0){
                                    for(u in 1:length(a_prendre)){
                                          ma_date_cur <- ma_date[a_prendre[[u]]]
                                          elevation_cur <- elevation[a_prendre[[u]]]
                                          snr_cur <- snr[a_prendre[[u]]]
                                          azimuth_cur <- azimuth[a_prendre[[u]]]

                                          date_base[[incr]] <- ma_date_cur
                                          nom_sat_base[incr] <- nom_sat
                                          elevation_base[[incr]] <- elevation_cur
                                          azimuth_base[[incr]] <- azimuth_cur
                                          SNR_base[[incr]] <- snr_cur
                                          incr <- incr + 1
                                    }
                              }
                        }
                  }
            }
      }
}

f <- NaN
f_log0 <- NaN
f_log1 <- NaN
incr_final <- 1
for(u in 1:length(SNR_base)){
      SNR_base_cur <- SNR_base[[u]]
      date_base_cur <- date_base[[u]]
      elevation_base_cur <- elevation_base[[u]]
      index_log0 <- seq(1,length(SNR_log0),1)[nom_sat_log0==nom_sat_base[[u]]]
      SNR_log0_tot <- NaN
      date_log0_tot <- NaN
      elevation_log0_tot <- NaN
      azimuth_log0_tot <- NaN
      if(length(index_log0)>=1){
            for(uu in 1:length(index_log0)){
                SNR_log0_tot <- c(SNR_log0_tot,SNR_log0[[index_log0[uu]]])
                date_log0_tot <- c(date_log0_tot,date_log0[[index_log0[uu]]])
                elevation_log0_tot <- c(elevation_log0_tot,elevation_log0[[index_log0[uu]]])
                azimuth_log0_tot <- c(azimuth_log0_tot,azimuth_log0[[index_log0[uu]]])
            }
            SNR_log0_cur <- SNR_log0_tot[2:length(SNR_log0_tot)]
            date_log0_cur <- date_log0_tot[2:length(date_log0_tot)]
            elevation_log0_cur <- elevation_log0_tot[2:length(elevation_log0_tot)]
            azimuth_log0_cur <- azimuth_log0_tot[2:length(azimuth_log0_tot)]

      }


      index_log1 <- seq(1,length(SNR_log1),1)[nom_sat_log1==nom_sat_base[[u]]]
      SNR_log1_tot <- NaN
      date_log1_tot <- NaN
      elevation_log1_tot <- NaN
      azimuth_log1_tot <- NaN
      if(length(index_log1)>=1){
            for(uu in 1:length(index_log1)){
                SNR_log1_tot <- c(SNR_log1_tot,SNR_log1[[index_log1[uu]]])
                date_log1_tot <- c(date_log1_tot,date_log1[[index_log1[uu]]])
                elevation_log1_tot <- c(elevation_log1_tot,elevation_log1[[index_log1[uu]]])
                azimuth_log1_tot <- c(azimuth_log1_tot,azimuth_log1[[index_log1[uu]]])
            }
            SNR_log1_cur <- SNR_log1_tot[2:length(SNR_log1_tot)]
            date_log1_cur <- date_log1_tot[2:length(date_log1_tot)]
            elevation_log1_cur <- elevation_log1_tot[2:length(elevation_log1_tot)]
            azimuth_log1_cur <- azimuth_log1_tot[2:length(azimuth_log1_tot)]

      }
      
      date_log0_cur <<- date_log0_cur[!is.na(SNR_log0_cur)]
      date_log1_cur <<- date_log1_cur[!is.na(SNR_log1_cur)]
      date_base_cur <<- date_base_cur[!is.na(SNR_base_cur)]
      SNR_log0_cur <<- SNR_log0_cur[!is.na(SNR_log0_cur)]
      SNR_log1_cur <<- SNR_log1_cur[!is.na(SNR_log1_cur)]
      SNR_base_cur <<- SNR_base_cur[!is.na(SNR_base_cur)]

      if(length(date_log0_cur)==0){
        date_log0_cur <- NaN
      }
      if(length(date_log1_cur)==0){
        date_log1_cur <- NaN
      }
      if(length(date_base_cur)==0){
        date_base_cur <- NaN
      }

 
      if(nom_x=='date'){
            date_min <- max(min(date_log0_cur),min(date_log1_cur),min(date_base_cur),na.rm=TRUE)
            date_max <- min(max(date_log0_cur),max(date_log1_cur),max(date_base_cur),na.rm=TRUE)
            ech_max <- min(mean(diff(date_log0_cur)),mean(diff(date_log1_cur)),mean(diff(date_base_cur)),na.rm=TRUE)
            if(date_max < date_min){
                date_min <- date_max
            }
            x <- seq(date_min,date_max,ech_max)

            SNR_log0_cur_reech <- NaN
            if(length(SNR_log0_cur)!=0){
                  SNR_log0_cur_reech <- approx(x=date_log0_cur,y=SNR_log0_cur,x)$y
            }

            SNR_log1_cur_reech <- NaN
            cor1 <- NaN
            if(length(SNR_log1_cur)!=0){
                  SNR_log1_cur_reech <- approx(x=date_log1_cur,y=SNR_log1_cur,xout=x)$y
            }

            SNR_base_cur_reech <- NaN
            cor2 <- NaN
            if(length(SNR_base_cur)!=0){
                  SNR_base_cur_reech <- approx(date_base_cur,SNR_base_cur,xout=x)$y
            }

      }

      if(nom_x=='sin'){
            x_min <- max(min(sin(elevation_log0_cur)),min(sin(elevation_log1_cur)),min(sin(elevation_base_cur)),na.rm=TRUE)
            x_max <- min(max(sin(elevation_log0_cur)),max(sin(elevation_log1_cur)),max(sin(elevation_base_cur)),na.rm=TRUE)
            ech_max <- min(mean(abs(diff(sin(elevation_log0_cur)))),mean(abs(diff(sin(elevation_log1_cur)))),mean(abs(diff(sin(elevation_base_cur)))),na.rm=TRUE)
            if(x_max < x_min){
                x_min <- x_max
            }
            x <- seq(x_min,x_max,ech_max)
            
            if(length(SNR_log0_cur)!=0){
                  SNR_log0_cur_reech <- approx(sin(elevation_log0_cur),y=SNR_log0_cur,xout=x)$y
            }
 
            SNR_log1_cur_reech <- NaN
            cor1 <- NaN
            if(length(SNR_log1_cur)!=0){
                  SNR_log1_cur_reech <- approx(sin(elevation_log1_cur),y=SNR_log1_cur,xout=x)$y
            }

            SNR_base_cur_reech <- NaN
            cor2 <- NaN
            if(length(SNR_base_cur)!=0){
                  SNR_base_cur_reech <- approx(sin(elevation_base_cur),SNR_base_cur,xout=x)$y
            }
      }
      
      y_min <- min(SNR_log0_cur_reech,SNR_log1_cur_reech,SNR_base_cur_reech,na.rm=TRUE)-30
      y_max <- max(SNR_log0_cur_reech,SNR_log1_cur_reech,SNR_base_cur_reech,na.rm=TRUE)

      if(length(x)>1){
            

      #      f_log0[incr_final] <- calcul_lsp(x,SNR_log0_cur_reech,2*h_min/longueur_onde,2*h_max/longueur_onde)*longueur_onde/2
#            f_log1[incr_final] <- calcul_lsp(x,SNR_log1_cur_reech,2*h_min/longueur_onde,2*h_max/longueur_onde)*longueur_onde/2
#            f[incr_final] <- calcul_lsp(x,SNR_log0_cur_reech-SNR_log1_cur_reech,2*h_min/longueur_onde,2*h_max/longueur_onde)*longueur_onde/2
            a <- scan('log0.pos',what=character(0),sep=';')
            a_mat <- matrix(a,ncol=15,byrow=TRUE)
            he <- as.numeric(a_mat[,5])
            lat <- as.numeric(a_mat[,3])
            long <- as.numeric(a_mat[,2])
            date_drone <- as.POSIXct(paste(a_mat[,1],a_mat[,2],sep=' '),origin='1970-01-01', tz='UTC')
             
            x11()
            plot(as.POSIXct(x,origin='1970-01-01',tz='UTC'),SNR_base_cur_reech,col='black',type='l',ylim=c(y_min,y_max),main=paste(nom_sat_base[[u]]),xlab='Date (2016-01-07 11h)',ylab='SNR (dBHz)',xlim=c(as.double(as.POSIXct('2016-01-07 11:00:00',origin='1970-01-01',tz='UTC')),as.double(as.POSIXct('2016-01-07 11:40:00',origin='1970-01-01',tz='UTC'))))
            flag <- TRUE
            if(length(SNR_log0_cur_reech)==length(SNR_base_cur_reech) & length(SNR_base_cur_reech)>longueur_minimum_serie){
                  lines(as.POSIXct(x,origin='1970-01-01',tz='UTC'),SNR_log0_cur_reech,col='red')
                  #text(mean(as.POSIXct(x,origin='1970-01-01',tz='UTC')),mean(SNR_log0_cur_reech),cor(SNR_base_cur_reech,SNR_log0_cur_reech),col='green',cex=2)
                  flag <- FALSE
            }
            if(length(SNR_log1_cur_reech)==length(SNR_base_cur_reech) & length(SNR_base_cur_reech)>longueur_minimum_serie){
                  lines(as.POSIXct(x,origin='1970-01-01',tz='UTC'),SNR_log1_cur_reech,col='blue')
                 # text(mean(as.POSIXct(x,origin='1970-01-01',tz='UTC')),mean(SNR_log1_cur_reech),cor(SNR_base_cur_reech,SNR_log1_cur_reech),col='green',cex=2)
                  if(flag==FALSE){
                        #text(mean(as.POSIXct(x,origin='1970-01-01',tz='UTC')),mean(SNR_base_cur_reech),cor(SNR_log0_cur_reech,SNR_log1_cur_reech),col='blue',cex=2)
                        y_max <- y_max*0.3
                        y_min <- y_min*1.8
                        lines(date_drone,(lat-min(lat))*(y_max-y_min)/(max(lat)-min(lat))+y_min+5,col='purple',type='l')
                        
                        lines(date_drone,(he-min(he))*(y_max-y_min)/(max(he)-min(he))+y_min+5,col='black',type='l')
                        
                        x11()
                        plot(as.POSIXct(x,origin='1970-01-01',tz='UTC'),SNR_log1_cur_reech/(SNR_log1_cur_reech+SNR_log0_cur_reech),,xlab='Date (2016-01-07 11h)',ylab='SNR',col='blue',type='l',main=paste(nom_sat_base[[u]]),ylim=c(-0.35,1.2),xlim=c(as.double(as.POSIXct('2016-01-07 11:00:00',origin='1970-01-01',tz='UTC')),as.double(as.POSIXct('2016-01-07 11:40:00',origin='1970-01-01',tz='UTC'))))
                        lines(as.POSIXct(x,origin='1970-01-01',tz='UTC'),SNR_log1_cur_reech/(SNR_log0_cur_reech),col='red',type='l')
                        y_max <- 1
                        y_min <- 0.3
                        lines(date_drone,(lat-min(lat))*(y_max-y_min)/(max(lat)-min(lat))-0.3,col='purple',type='l')
                        lines(date_drone,(he-min(he))*(y_max-y_min)/(max(he)-min(he))-0.3,col='black',type='l')
                  }
            }
            else{
                  if(flag==TRUE){
                      dev.off()
                  }
            }
      }
      
      incr_final <- incr_final + 1
}


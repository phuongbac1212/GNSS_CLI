calcul_coord <- function(mon_data,sp3_p,sp3_c,sp3_s){

      # Lecture des fichier sp3
      fichier <- file(paste(path_input,sp3_p,sep=''), open = "r")
      lines_p <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",readLines(fichier,warn=FALSE))),' ')
      nb_info <- as.numeric(lines_p[[3]][2])
      close(fichier)
      lines_date <- lines_p[substr(sapply(lines_p,'[[',1),1,1)=='*']
      annee <- sapply(lines_date,'[[',2)
      mois <- sapply(lines_date,'[[',3)
      jour <- sapply(lines_date,'[[',4)
      heure <- sapply(lines_date,'[[',5)
      minute <- sapply(lines_date,'[[',6)
      seconde <- sapply(lines_date,'[[',7)
      
      date_i <- as.double(strptime(paste(annee,'-',mois,'-',jour,' ',heure,':',minute,':',seconde,sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'))-17
      date_sp3 <- rep(date_i[1],nb_info)
      for(u in 2:length(date_i)){
            date_sp3 <- c(date_sp3,rep(date_i[u],nb_info))
      }
      
      fichier <- file(paste(path_input,sp3_c,sep=''), open = "r")
      lines_c <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",readLines(fichier))),' ')
      nb_info <- as.numeric(lines_c[[3]][2])
      close(fichier)
      lines_date <- lines_c[substr(sapply(lines_c,'[[',1),1,1)=='*']
      annee <- sapply(lines_date,'[[',2)
      mois <- sapply(lines_date,'[[',3)
      jour <- sapply(lines_date,'[[',4)
      heure <- sapply(lines_date,'[[',5)
      minute <- sapply(lines_date,'[[',6)
      seconde <- sapply(lines_date,'[[',7)
      date_i <- as.double(strptime(paste(annee,'-',mois,'-',jour,' ',heure,':',minute,':',seconde,sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'))-17
      for(u in 1:length(date_i)){
            date_sp3 <- c(date_sp3,rep(date_i[u],nb_info))
      }
      
      fichier <- file(paste(path_input,sp3_s,sep=''), open = "r")
      lines_s <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",readLines(fichier))),' ')
      nb_info <- as.numeric(lines_s[[3]][2])
      close(fichier)      
      lines_date <- lines_s[substr(sapply(lines_s,'[[',1),1,1)=='*']
      annee <- sapply(lines_date,'[[',2)
      mois <- sapply(lines_date,'[[',3)
      jour <- sapply(lines_date,'[[',4)
      heure <- sapply(lines_date,'[[',5)
      minute <- sapply(lines_date,'[[',6)
      seconde <- sapply(lines_date,'[[',7)
      date_i <- as.double(strptime(paste(annee,'-',mois,'-',jour,' ',heure,':',minute,':',seconde,sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'))-17
      for(u in 1:length(date_i)){
            date_sp3 <- c(date_sp3,rep(date_i[u],nb_info))
      }
      liness  <- c(lines_p,lines_c,lines_s)
      lines_data <- liness[substr(sapply(liness,'[[',1),1,1)=='P']
      lines_date <- date_sp3

      for(num_sat in 1:length(mon_data)){
            taille <- length(mon_data[[num_sat]])
            nom_sat <- mon_data[[num_sat]][1]
            x_sat <- sapply(lines_data,'[[',2)
            y_sat <- sapply(lines_data,'[[',3)
            z_sat <- sapply(lines_data,'[[',4)
            nom_sat_sp3 <- substr(sapply(lines_data,'[[',1),2,4)
            x_sat <- as.numeric(x_sat[nom_sat_sp3 == nom_sat])
            y_sat <- as.numeric(y_sat[nom_sat_sp3 == nom_sat])
            z_sat <- as.numeric(z_sat[nom_sat_sp3 == nom_sat])
            date_sat <- as.numeric(lines_date[nom_sat_sp3 == nom_sat])
            datee <- as.numeric(mon_data[[num_sat]][[2]])
            offsett <- 7200
            if(length(datee)!=0){
                  if(length(x_sat)!=0){
                        if(length(x_sat)!=0 & length(date_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)])!=0){ #& !is.na(date_sat < (max(datee,rm.na=TRUE)+offsett))
                              mon_data[[num_sat]][[taille+1]] <- spline(date_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)],x_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)],xout=datee)$y*1000
                              mon_data[[num_sat]][[taille+2]] <- spline(date_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)],y_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)],xout=datee)$y*1000
                              mon_data[[num_sat]][[taille+3]] <- spline(date_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)],z_sat[date_sat>(min(datee,rm.na=TRUE)-offsett) & date_sat < (max(datee,rm.na=TRUE)+offsett)],xout=datee)$y*1000            
                        }   
                        else{
                              mon_data[[num_sat]][[taille+1]] <- NaN
                              mon_data[[num_sat]][[taille+2]] <- NaN
                              mon_data[[num_sat]][[taille+3]] <- NaN
                        }
                  }
            }

      }
      
      return(mon_data)
}
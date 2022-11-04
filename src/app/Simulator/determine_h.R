determine_h <- function(autonome){
     library(MASS,quietly=TRUE)
      if(autonome){
            source('constants.r') 
            source('lecture_directeur.r') 
      }
      
      mon_data <-matrix(scan(paste(path_output,'sortie_f.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=14,byrow=TRUE)
      
      # on applique les filtres en elevations
      mon_data <- mon_data[as.numeric(mon_data[,8])>elevation_min & as.numeric(mon_data[,8])<elevation_max,]
      mon_data <- mon_data[!is.na(mon_data[,1]),]

      #on applique les filtres en azimuth
      if(dim(mon_data)[1]!=0){
            boolean <- seq(1,dim(mon_data)[1],1)*0  
            seq_a_prendre <- seq(1,dim(mon_data)[1],1)    
            azi <- as.numeric(mon_data[,9])
            for(u in seq(1,length(masque_azimuth),2)){
                  boolean[(as.numeric(mon_data[,9])>=masque_azimuth[u] & as.numeric(mon_data[,9]) <= masque_azimuth[u+1])] <- boolean[(as.numeric(mon_data[,9])>=masque_azimuth[u] & as.numeric(mon_data[,9]) <= masque_azimuth[u+1])]+1
            }
            boolean[boolean!=0] <- 1
            seq_a_prendre <- seq_a_prendre[as.logical(boolean)]
            mon_data <- mon_data[seq_a_prendre,]
             
             # on prend que les satellites qui nous interessent           
            nom_sat <- mon_data[,7]
            ma_date <- as.double(as.POSIXct(paste(mon_data[,1],'-',mon_data[,2],'-',mon_data[,3],' ',mon_data[,4],':',mon_data[,5],':',mon_data[,6],sep=''),format="%Y-%m-%d %H:%M:%S",tz='UTC'),tz='UTC')[nom_sat %in% nom_sat_liste]
            f <- as.numeric(mon_data[,14])[nom_sat %in% nom_sat_liste]
            longueur_onde <- as.numeric(mon_data[,13])[nom_sat %in% nom_sat_liste]
            elevation <- as.numeric(mon_data[,8])[nom_sat %in% nom_sat_liste]
            elevation_point <- as.numeric(mon_data[,10])[nom_sat %in% nom_sat_liste]
            amplitude <- as.numeric(mon_data[,11])
            phase <- as.numeric(mon_data[,12])
      
            nb_equ_min <- 5   #nombre d'equations minimums dans la resolution par les moindres carres   
            filtrer_f_1_sigma <- FALSE
            filtrer_h_1_sigma <- FALSE
            
            ma_date_h <- NaN
            h <- NaN
            h_point <- NaN
            emq_h <- NaN
            emq_h_point <- NaN
            ma_date <- ma_date[!is.na(f)]
            longueur_onde <- longueur_onde[!is.na(f)]
            elevation <- elevation[!is.na(f)]
            elevation_point <- elevation_point[!is.na(f)]
            f <- f[!is.na(f)]
        
            if(filtrer_f_1_sigma){     
                  a <- mean(f)+sd(f)
                  b <- mean(f)-sd(f)
                  ma_date <- ma_date[f<a & f>b]
                  elevation <- elevation[f<a & f>b]
                  elevation_point <- elevation_point[f<a & f>b]
                  longueur_onde <- longueur_onde[f<a & f>b]
                  f <- f[f<a & f>b] 
            }
            if(length(f)!=0){
                  if(methode_a_utiliser=='Larson'){
                        ma_seq <- seq(1,length(ma_date),1)
                        ma_date_cur <- min(ma_date)+ecart_max
                        incr <- 1
                        while(ma_date_cur<=max(ma_date)){
                              a_prendre <- ma_seq[abs(ma_date-ma_date_cur)<=ecart_max]
                              ma_date_h[incr] <- ma_date_cur
                              h[incr] <- mean(f[a_prendre]*longueur_onde[a_prendre]/2)
                              ma_date_cur <- ma_date_cur+1/Fh
                              incr <- incr + 1
                        }
                        emq_h <- rep(NaN,length(h))
                        emq_h_point <- rep(NaN,length(h))
                        h_point <- rep(NaN,length(h))
                  }
            
                  if(methode_a_utiliser == 'LSM'){
            
                        liste_t <- seq(min(ma_date),max(ma_date),by=1/Fh) 
            
                        incr <- 1
                        for(t in liste_t){ # pour chaque instant t de la liste a calculer
                              a_prendre <- (abs(ma_date-t) <= ecart_max)
                              ecart <- abs(ma_date-t)
                              # on va resoudre K1*h + K2*h_point + K3 = 0 par les moindres carres            
                              K2 <- tan(elevation[a_prendre])/elevation_point[a_prendre]
                              K1 <- seq(from=1,to=1,length=length(K2))
                              K3 <- f[a_prendre]*longueur_onde[a_prendre]/2 
                              if(length(K3)>nb_equ_min){
                                    A <- mat.or.vec(length(K1),2)
                                    A[,1] <- K1
                                    A[,2] <- K2
                                    P <- diag(length(K1))
                                    mamat <- NA
                                    library(MASS)
                                    try(mamat <- solve(t(A)%*%P%*%A),silent=TRUE)
                                    if(length(mamat)>1){
                                          result <- mamat%*%(t(A)%*%P%*%K3)
                                          h[incr] <- result[1]
                                          h_point[incr] <- result[2]                           
                                          ma_date_h[incr] <- t
                                          v <- A%*%c(h[incr],h_point[incr])-K3
                                          emq_0 <- sqrt((t(v)%*%P%*%v)/(length(K1)-2))
                                          emq_h[incr] <- emq_0*sqrt(mamat[1,1])
                                          emq_h_point[incr] <- emq_0*sqrt(mamat[2,2]) 
                                          incr <- incr + 1                           
                                    }
                              }
                              
                        }
                        
                 } # fin de methode LSM   
            
                  ma_date_h <- ma_date_h[!is.na(h)]
                  emq_h <- emq_h[!is.na(h)]
                  emq_h_point <- emq_h[!is.na(h)]
                  h_point <- h_point[!is.na(h)]    
                  h <- h[!is.na(h)]
                  ma_date_h <- ma_date_h[h<h_max & h>h_min]
                  emq_h <- emq_h[h<h_max & h>h_min]
                  h_point <- h_point[h<h_max & h>h_min]
                  emq_h_point <- emq_h_point[h<h_max & h>h_min]
                  h <- h[h<h_max & h>h_min] 
                  if(filtrer_h_1_sigma){     
                        a <- mean(h)+sd(h)
                        b <- mean(h)-sd(h)
                        ma_date_h <- ma_date_h[h<a & h>b]
                        emq_h <- emq_h[h<a & h>b]
                        emq_h_point <- emq_h_point[h<a & h>b]
                        h_point <- h_point[h<a & h>b]
                        h <- h[h<a & h>b] 
                  }
            }
            else{
                  print('Pas assez de f dans le fichier sortie_f.txt pour pouvoir calculer h',quote=FALSE)
                  h <- rep(NaN,length(ma_date_h))
                  emq_h <- rep(NaN,length(ma_date_h))
                  h_point <- rep(NaN,length(ma_date_h))
                  
            } 
            library(lubridate)      
            date_final <- as.POSIXct(ma_date_h,origin='1970-01-01',tz='UTC')
            col1 <- c('Year',format(as.numeric(year(as.POSIXct(ma_date_h,origin='1970-01-01',tz='UTC'))),digits=1))
            col2 <- c('Month',month(date_final))
            col3 <- c('Day',day(date_final))
            col4 <- c('Hour',hour(date_final))
            col5 <- c('Minute',minute(date_final))
            col6 <- c('Seconde',format(second(date_final),digits=1))
            col7 <- c('h',format(h,digits=8))
            col8 <- c('emq_h',format(emq_h,digits=8))
            col9 <- c('h_point',format(h_point,digits=8))
            col10 <- c('emq_h_point',format(emq_h_point,digits=8))
            fichier_sortie <- paste(path_output,'sortie_h.txt',sep='')
            write.matrix(matrix(c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10),ncol=10), file = fichier_sortie,sep='\t')  
            #plot(as.POSIXct(ma_date_h,origin='1970-01-01'),h,type='l')
#            print(paste('Moyenne =',mean(h),'m'))
#            print(paste('Ecart-type =',sd(h),'m'))     
#            print(paste('Mediane = ',median(h),'m'))
            
      }
}









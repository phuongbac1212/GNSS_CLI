calcul_azimuth_elevation_et_enregistre <- function(nom_sat_liste,mon_dataa){

      for(num_sat in 1:length(nom_sat_liste)){
            longueur <- length(mon_dataa[[num_sat]])
            if(length(mon_dataa[[num_sat]][[longueur]])>10){
                  date_cur <- mon_dataa[[num_sat]][[2]]
                  pos_R_geo_cur <- cbind(pos_R_geo[,1]*180/pi,pos_R_geo[,2]*180/pi,pos_R_geo[,3])
                  if(!recepteur_statique){
                        long_R <- spline(x=date_rec,y=pos_R_geo[,1],xout=date_cur,rule=2)$y
                        lat_R <- spline(x=date_rec,y=pos_R_geo[,2],xout=date_cur,rule=2)$y
                        he_R <- spline(x=date_rec,y=pos_R_geo[,3],xout=date_cur,rule=2)$y
                        pos_R_geo_cur <- cbind(long_R,lat_R,he_R)
                  }
                  pos_R_xyz_cur <- geo2pos(R_0_WGS84,R_p_WGS84,pos_R_geo_cur) 
                  pos_R_0_xyz <- geo2pos(R_0_WGS84,R_p_WGS84,(pos_R_geo_cur+c(rep(0,dim(pos_R_geo_cur)[1]),rep(0,dim(pos_R_geo_cur)[1]),rep(-5,dim(pos_R_geo_cur)[1]))))
                  vec_rec_0 <- rbind(pos_R_0_xyz - pos_R_xyz_cur)
                  X_sat <- mon_dataa[[num_sat]][[longueur-2]]
                  Y_sat <- mon_dataa[[num_sat]][[longueur-1]]
                  Z_sat <- mon_dataa[[num_sat]][[longueur]]
                  pos_sat_xyz <- cbind(X_sat,Y_sat,Z_sat)
                  coord_sphere_sat <- cart2sph(pos_sat_xyz)*180/pi
                  coord_sphere_rec <- cart2sph(pos_R_xyz_cur)*180/pi
                  coord_sphere_rec <<- coord_sphere_rec

                  if(recepteur_statique){
                        coord_sphere_rec <- matrix(coord_sphere_rec,ncol=3)
                  }
                  azimuth <- (bearing(cbind(coord_sphere_rec[,1],coord_sphere_rec[,2]),cbind(coord_sphere_sat[,1],coord_sphere_sat[,2])))*pi/180 
                  mon_dataa[[num_sat]][[longueur+1]] <- azimuth
                  
                  vec_rec_sat <- cbind(X_sat-pos_R_xyz_cur[,1],Y_sat-pos_R_xyz_cur[,2],Z_sat-pos_R_xyz_cur[,3])
                  scalaire <- vec_rec_sat[,1]*vec_rec_0[,1] + vec_rec_sat[,2]*vec_rec_0[,2] + vec_rec_sat[,3]*vec_rec_0[,3]
      
                  cos_elevation_1 <- scalaire/sqrt(vec_rec_0[,3]^2+vec_rec_0[,2]^2+vec_rec_0[,1]^2)
                  cos_elevation_2 <- cos_elevation_1/sqrt(vec_rec_sat[,3]^2+vec_rec_sat[,2]^2+vec_rec_sat[,1]^2)
                  acos_elevation <- acos(cos_elevation_2)
                  elevation <- acos_elevation - pi/2
                  elevation[elevation>=90*pi/180] <- 180*pi/180-elevation[elevation>=90*pi/180]
                  mon_dataa[[num_sat]][[longueur+2]] <-  elevation

                  # on calcule la distance sat <--> rec

                  distance <- sqrt(apply(vec_rec_sat^2,1,sum))
                  mon_dataa[[num_sat]][[longueur+3]] <- distance
                  nom_sat_cur <- gsub(' ','0',mon_dataa[[num_sat]][1])
                  fichier_sortie <- paste(path_extraction,'extraction_',rinex_cur,'_',nom_sat_cur,'.tmp',sep='')
                  nb_SNR <- length(mon_dataa[[num_sat]])-9
                  SNR <- c(mon_dataa[[num_sat]][[length(mon_dataa[[num_sat]])-6]][1],mon_dataa[[num_sat]][[3]])
                  if(nb_SNR >=2){
                        for(u in 2:nb_SNR){
                              SNR <- c(SNR,mon_dataa[[num_sat]][[length(mon_dataa[[num_sat]])-6]][u],mon_dataa[[num_sat]][[u+2]])
                        }
                  }
  
                  SNR[SNR==''] <- '-----'
                  SNR[is.na(SNR)] <- '-----'
                  SNR[SNR<=1 & SNR>=0] <- '-----'
                  library(lubridate)
                  date_cur <- as.POSIXct(as.numeric(mon_dataa[[num_sat]][[2]]),tz='UTC',origin='1970-01-01')
                  col1 <- c('y',format(year(as.POSIXct(as.numeric(mon_dataa[[num_sat]][[2]]),tz='UTC',origin='1970-01-01')),digits=1))
                  col2 <- c('m',month(date_cur))
                  col3 <- c('d',day(date_cur))
                  col4 <- c('h',hour(date_cur))
                  col5 <- c('m',minute(date_cur))
                  col6 <- c('s',format(second(date_cur),digits=1))
                  col7 <- c('PRN',rep(nom_sat_cur,length(date_cur)))
                  col8 <- c('elevation(rad)',format(elevation,digits=8))
                  col9 <- c('azimuth(rad)',format(azimuth,digits=8))
                  col10 <- c('distance(m)',format(distance,digits=10))                
                  write.matrix(matrix(c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,format(SNR,digits=3)),ncol=(10+nb_SNR)), file = fichier_sortie,sep='\t')
            }
      }
}
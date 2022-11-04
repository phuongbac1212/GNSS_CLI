##################################################################
#
# janvier 2016 - Nicolas ROUSSEL
#
# packages necessaires : pracma, geosphere, sp, lubridate, MASS, lomb, tcltk  --> faire install.packages('nom_du_package')
#
##################################################################


rm(list=ls())  # vide le cache
graphics.off()  # ferme toutes les figures deja ouvertes
options(warn=-1) # Naffiche pas les warnings
print(' ',quote=FALSE)
print(' ',quote=FALSE)

######## Definition des constantes ########
source('constants.r') 
      
# On lit le fichier directeur et on sort toutes les infos
source('lecture_directeur.r') 




############################
# Si l'utilisateur a choisit de lancer le programme EXTRACTION
############################


if(faire_EXTRACTION){
      print('-------------------------------------------------------------------------------------------------------------------------------------------------------------',quote=FALSE)
      print('************************************************************** I.EXTRACTION **************************************************************',quote=FALSE)
      print('-------------------------------------------------------------------------------------------------------------------------------------------------------------',quote=FALSE)
      T_debut_extraction <- Sys.time()
      
      
      # Importation des fonctions necessaires 
      library(pracma,quietly=TRUE)
      library(geosphere,quietly=TRUE)
      library(lubridate,quietly=TRUE)
      library(MASS,quietly=TRUE)
      source('determine_sp3_correspondant.r')
      source('jul_2_gps.r')
      source('extraction_snr.r')
      source('telecharger_sp3.r')
      source('calcul_coord.r')
      source('lecture_donnees.r')
      source('calcul_azimuth_elevation_et_enregistre.r')     
      
      liste_rinex <- c(list.files(path = path_input, pattern = ("o$")),list.files(path = path_input, pattern = ("obs$")))  
      
      if(length(liste_rinex)==0){
            print('Il n\'y a aucun fichier rinex a extraire dans le dossier 0.input',quote=FALSE); 
      }
      for(num_rinex in 1:length(liste_rinex)){   # pour chaque fichier d'observation present dans le dossier 0.input
        rinex_cur <-liste_rinex[num_rinex]
        fichier <- file(paste(path_input,rinex_cur,sep=''), open = "r")
        print('',quote=FALSE)
        print(paste('**** Extraction en cours du fichier : ',rinex_cur,' ****',sep=''),quote=FALSE)
        while (length(line <- scan(fichier, what = character(0), sep = "", nlines = 1, quiet = TRUE)) > 0) { # tant qu'on n'arrive pas a la fin
          long <- length(line)
          if(line[long]=='OBS' & line[long-1]=='FIRST'){
            annee_debut <- line[1]
            mois_debut <- line[2]
            jour_debut <- as.numeric(line[3])
            heure_debut <- line[4]
            minute_debut <- line[5]
            seconde_debut <- line[6] 
          }
          if(line[long]=='OBS' & line[long-1]=='LAST'){
            annee_fin <- line[1]
            mois_fin <- line[2]
            jour_fin <- as.numeric(line[3])
            heure_fin <- line[4]
            minute_fin <- line[5]
            seconde_fin <- line[6]
            break
          }
        }
        
  
        extraction(rinex_cur,annee_debut,mois_debut,jour_debut,annee_fin,mois_fin,jour_fin)
        print('',quote=FALSE)
        print('',quote=FALSE)
        print(paste('Lextraction est finie'),quote=FALSE) 
        T_fin_extraction <- Sys.time() 
        print(paste('Temps d\'execution de l\'extraction :',difftime(T_fin_extraction,T_debut_extraction)),quote='FALSE')
        print('',quote=FALSE)
        print('',quote=FALSE)
        print('',quote=FALSE)
      }
      closeAllConnections()
}














############################
# Si lutilisateur a choisi de lancer le programme ANALYSE
############################

if(faire_ANALYSE){
      library(MASS,quietly=TRUE)
      print('-------------------------------------------------------------------------------------------------------------------------------------------------------------',quote=FALSE)
      print('************************************************************** II.ANALYSE **************************************************************',quote=FALSE)
      print('-------------------------------------------------------------------------------------------------------------------------------------------------------------',quote=FALSE)
      
      
      T_debut_f <- Sys.time()
      option <- 0
      
      source('determine_f.r')
      if(system_exploitation == "windows"){
            system('cmd /c del_output.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
      }
      if(system_exploitation == "linux"){
            system('rm ../2.output/*.txt')
            system('rm ../2.output/*.tmp')
      }
      # on verifie que le dossier a traiter nest pas vide
      liste_extraction_complet <- list.files(path = path_extraction, pattern = (".tmp"))# on va traiter l'ensemble des fichiers presents dans le dossier 1.extraction
      if(length(liste_extraction_complet)==0){
            print('Il n\'y a aucun fichier \'extraction\' a traiter dans le dossier 1.extraction',quote=FALSE); 
      }
      
      # pour chaque fichier dans le dossier dextraction, on determine f, A, phi
      for(fichier in liste_extraction_complet){
            determine_f(fichier,option)                                             
      }
      
      source('concat_fichiers_f.r')   
      
      T_fin_f <-Sys.time() 
      closeAllConnections()
      
      # Calcul de h  
      print(' ',quote=FALSE)
      print('Determination de h',quote=FALSE)
      T_debut_h <- Sys.time()
      option <- 0
      source('determine_h.r')
      liste_f_complet <- list.files(path = path_output, pattern = (".txt"))# on va traiter le fichier present dans le dossier 2.output
      if(length(liste_f_complet)==0 | !('sortie_f.txt' %in% liste_f_complet)){
            print('Il n\'y a aucun fichier \'sortie_f.txt\' a traiter dans le dossier 2.output',quote=FALSE); 
      }
      else{
            source('determine_h.r');determine_h(option)
      }    
      closeAllConnections()
      T_fin_h <- Sys.time()
      print('',quote=FALSE)
      print('',quote=FALSE)
      print('Lanalyse est finie',quote=FALSE)
     # print(paste('Temps dexecution de lanalyse:',difftime(T_fin_h,T_debut_f)),quote='FALSE')
      print('',quote=FALSE)
      print('',quote=FALSE)
      print('',quote=FALSE)
}












############################
# Cette partie est uniquement experimentale --> non validee !!!!!
############################

faire_analyse_resultats <- FALSE
if(faire_analyse_resultats){

      moyenne_glissante <- TRUE
      delta_t <- 28800 # en s
      T_debut_analyse_resultats <- Sys.time()
      source('clustering.r')
      angle_coupure <- 30*pi/180

      liste_f_complet <- list.files(path = path_output, pattern = (".txt"))# on va traiter le fichier present dans le dossier 2.output
      if(length(liste_f_complet)==0 | !('sortie_f.txt' %in% liste_f_complet)){
            print('Il n\'y a aucun fichier \'sortie_f.txt\' a traiter dans le dossier 2.output',quote=FALSE); 
      }
      
      if(length(liste_f_complet)==0 | !('sortie_h.txt' %in% liste_f_complet)){
            print('Il n\'y a aucun fichier \'sortie_h.txt\' a traiter dans le dossier 2.output',quote=FALSE); 
      }
      mon_data_f <-matrix(scan(paste(path_output,'sortie_f.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=14,byrow=TRUE)
      mon_data_h <-matrix(scan(paste(path_output,'sortie_h.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=10,byrow=TRUE)    
      
      # on prend que les satellites qui nous interessent           
      nom_sat_f <- mon_data_f[,7]
      ma_date_h <- as.double(as.POSIXct(paste(mon_data_h[,1],'-',mon_data_h[,2],'-',mon_data_h[,3],' ',mon_data_h[,4],':',mon_data_h[,5],':',mon_data_h[,6],sep=''),format="%Y-%m-%d %H:%M:%S",tz='UTC'),tz='UTC')
      h <- as.numeric(mon_data_h[,7])
      hh <- 0
          for(u in 1:length(h)){
                                          hh[u] <- mean(h[ma_date_h<(ma_date_h[u]+delta_t) & ma_date_h>(ma_date_h[u]-delta_t)],na.rm=TRUE)      
                                    }  
                                    
      plot(as.POSIXct(paste(mon_data_h[,1],'-',mon_data_h[,2],'-',mon_data_h[,3],' ',mon_data_h[,4],':',mon_data_h[,5],':',mon_data_h[,6],sep=''),format="%Y-%m-%d %H:%M:%S",tz='UTC'),hh)

#      azd
#      plot(as.POSIXct(paste(mon_data_f[,1],'-',mon_data_f[,2],'-',mon_data_f[,3],' ',mon_data_f[,4],':',mon_data_f[,5],':',mon_data_f[,6],sep=''),format="%Y-%m-%d %H:%M:%S",tz='UTC')[!is.na(mon_data_f[,14])],as.numeric(mon_data_f[,14][!is.na(mon_data_f[,14])])*as.numeric(mon_data_f[,13][!is.na(mon_data_f[,14])])/2)
##      plot(as.POSIXct(paste(mon_data_f[,1],'-',mon_data_f[,2],'-',mon_data_f[,3],' ',mon_data_f[,4],':',mon_data_f[,5],':',mon_data_f[,6],sep=''),format="%Y-%m-%d %H:%M:%S",tz='UTC'),mon_data_f[,11])
#       zef
      if(methode_a_utiliser=='Larson'){
            serie_concat <- list()
            serie_concat[[1]] <- NaN
            serie_concat[[2]] <- NaN
            serie_concat[[3]] <- NaN
            ma_date_concat <- list()
            ma_date_concat[[1]] <- NaN
            ma_date_concat[[2]] <- NaN
            ma_date_concat[[3]] <- NaN
            serie_finale <- list()
            serie_finale[[1]] <- NaN
            serie_finale[[2]] <- NaN
            serie_finale[[3]] <- NaN
            ma_date_finale <- list()
            ma_date_finale[[1]] <- NaN
            ma_date_finale[[2]] <- NaN
            ma_date_finale[[3]] <- NaN
            for(sat_a_utiliser in nom_sat_liste){
                  amplitude <- as.numeric(mon_data_f[,11])[nom_sat_f == sat_a_utiliser]
                  phase <- as.numeric(mon_data_f[,12])[nom_sat_f == sat_a_utiliser]
                  ma_date <- as.double(as.POSIXct(paste(mon_data_f[,1],'-',mon_data_f[,2],'-',mon_data_f[,3],' ',mon_data_f[,4],':',mon_data_f[,5],':',mon_data_f[,6],sep=''),format="%Y-%m-%d %H:%M:%S",tz='UTC'),tz='UTC')[nom_sat_f == sat_a_utiliser]          
                  elevation <- as.numeric(mon_data_f[,8])[nom_sat_f == sat_a_utiliser]
                  azimuth <- as.numeric(mon_data_f[,9])[nom_sat_f == sat_a_utiliser]
                  f <- as.numeric(mon_data_f[,14])[nom_sat_f == sat_a_utiliser]
                  longueur_onde <- as.numeric(mon_data_f[,13])[nom_sat_f == sat_a_utiliser]
                  elevation_point <- as.numeric(mon_data_f[,10])[nom_sat_f == sat_a_utiliser]
                  h <- f*longueur_onde/2


                  amplitude <- amplitude[!is.na(f)]
                  ma_date <- ma_date[!is.na(f)]
                  elevation <- elevation[!is.na(f)]
                  elevation_point <- elevation_point[!is.na(f)]
                  azimuth <- azimuth[!is.na(f)]
                  phase <- phase[!is.na(f)]
                  h <- h[!is.na(f)]
                  longueur_onde <- longueur_onde[!is.na(f)]

                  if(length(amplitude)>= 2){
                        a_prendre <- clustering(elevation_point,1,500,500,ma_date_sat,elevation_point,12)# on cherche les clusters (max 7)
                        for(serie_a_analyser in 1:3){
                              if(serie_a_analyser == 1){serie <- amplitude;titre <- 'A';signe <- 1}
                              else if(serie_a_analyser == 2){serie <- h;titre <- 'h';signe <- 1}
                              else if(serie_a_analyser == 3){serie <- phase;titre <- 'phi';signe <- -1}
                              for(num_cluster in 1:length(a_prendre)){
                                    serie_cur <- serie[a_prendre[[num_cluster]]]
                                    ma_date_cur <- ma_date[a_prendre[[num_cluster]]]
                                    elevation_cur <- elevation[a_prendre[[num_cluster]]]
                                    ma_date_cur <- ma_date_cur[!is.na(serie_cur)]
                                    elevation_cur <- elevation_cur[!is.na(serie_cur)]
                                    serie_cur <- serie_cur[!is.na(serie_cur)]
                                    if(length(serie_cur)>0){
                                          if(mean(elevation_cur,na.rm=TRUE)<angle_coupure){
                                                serie_cur <- signe*serie_cur
                                          }
                                          else{
                                                serie_cur <- -signe*serie_cur
                                          }

                                          serie_cur <- (serie_cur-min(serie_cur,na.rm=TRUE))*(1-0)/(max(serie_cur,na.rm=TRUE)-min(serie_cur,na.rm=TRUE))+0
                                          serie_concat[[serie_a_analyser]] <- c(serie_concat[[serie_a_analyser]],serie_cur)
                                          ma_date_concat[[serie_a_analyser]] <- c(ma_date_concat[[serie_a_analyser]],ma_date_cur)
                                         
                                    }
                              } # fin num_cluster
                              serie_concat[[serie_a_analyser]] <- serie_concat[[serie_a_analyser]][order(ma_date_concat[[serie_a_analyser]])]
                              ma_date_concat[[serie_a_analyser]] <- ma_date_concat[[serie_a_analyser]][order(ma_date_concat[[serie_a_analyser]])]
                              
                        }# fin serie_a_analyser  
                  }
            }
            for(serie_a_analyser in 1:3){
                  if(serie_a_analyser == 1){serie <- amplitude;titre <- 'A';signe <- 1}
                  else if(serie_a_analyser == 2){serie <- h;titre <- 'h';signe <- 1}
                  else if(serie_a_analyser == 3){serie <- phase;titre <- 'phi';signe <- -1}
                  ma_date_finale[[serie_a_analyser]] <- ma_date_concat[[serie_a_analyser]]
                              serie_finale[[serie_a_analyser]] <- serie_concat[[serie_a_analyser]]
                              
                              if(moyenne_glissante){
                                    for(u in 1:length(ma_date_concat[[serie_a_analyser]])){
                                          serie_finale[[serie_a_analyser]][u] <- mean(serie_concat[[serie_a_analyser]][ma_date_concat[[serie_a_analyser]]>(ma_date_concat[[serie_a_analyser]][u]-delta_t) & ma_date_concat[[serie_a_analyser]]<(ma_date_concat[[serie_a_analyser]][u]+delta_t)],na.rm=TRUE)      
                                          ma_date_finale[[serie_a_analyser]][u] <- ma_date_concat[[serie_a_analyser]][u]
                                    }      
                              }
                              else{
                                    serie_finale[[serie_a_analyser]] <- serie_concat[[serie_a_analyser]]
                                    ma_date_finale[[serie_a_analyser]] <- ma_date_concat[[serie_a_analyser]]
                              }
                              x11()
                              plot(as.POSIXct(ma_date_finale[[serie_a_analyser]],origin='1970-01-01'),serie_finale[[serie_a_analyser]],type='l',lwd=3,main=paste(titre,' = f(t)',sep=''),ylab=titre,xlab='Date')
                              datee <- as.POSIXct(ma_date_finale[[serie_a_analyser]],origin='1970-01-01')
                              
                              library(lubridate,quietly=TRUE)
                              library(MASS,quietly=TRUE)
                              col1 <- c('Year',year(datee))
                              col2 <- c('Month',month(datee))
                              col3 <- c('Day',day(datee))
                              col4 <- c('Hour',hour(datee))
                              col5 <- c('Minute',minute(datee))
                              col6 <- c('Second',second(datee))
                              col7 <- c(titre,serie_finale[[serie_a_analyser]])
                              write.matrix(matrix(c(col1,col2,col3,col4,col5,col6,col7),ncol=7), file = paste(path_output,'analyse_',titre,'.txt',sep=''),sep='\t')
            }
                 
      }# fin de if methode Larson 

      T_fin_analyse_resultats <- Sys.time() 
      print(paste('Temps d\'execution de l\'analyse des resultats :',difftime(T_fin_analyse_resultats,T_debut_analyse_resultats)))     
}









############################
# Si lutilisateur a choisi de lancer le simulateur
############################

if(faire_SIMULATEUR){
      # On importe les fonctions necessaires
      source('determine_sp3_correspondant.r')
      source('jul_2_gps.r')
      source('telecharger_sp3.r')
      library(MASS,quietly=TRUE)  


      print('-------------------------------------------------------------------------------------------------------------------------------------------------------------',quote=FALSE)
      print('************************************************************** IV.SIMULATION **************************************************************',quote=FALSE)
      print('-------------------------------------------------------------------------------------------------------------------------------------------------------------',quote=FALSE)
      
      
      T_debut_simu <- Sys.time()
      
      if(!recepteur_statique){ # si le recepteur est mobile
            liste_fichier_pos_complet <<- list.files(path = path_input, pattern = (".pos"))
            fichier_pos <- gsub(path_input,'',nom_fichier_pos)
            if(length(liste_fichier_pos_complet)==0 | !(fichier_pos %in% liste_fichier_pos_complet)){
                  print(paste('Erreur : il n\'y a aucun fichier ',fichier_pos, ' a traiter dans le dossier 0.input. Veuillez v?rifier.',sep=''),quote=FALSE); 
            }
            
            coord_rec <- scan(nom_fichier_pos,what=character(0),quiet=TRUE)
            coord_rec <- coord_rec[(which(coord_rec=='ratio') + 1) : length(coord_rec)]
            coord_rec <- matrix(coord_rec,ncol=15,byrow=TRUE)
            date_rec <-  as.double(as.POSIXct(paste(gsub('/','-',coord_rec[,1]),' ',coord_rec[,2],sep='')))
            lat_rec <- as.numeric(coord_rec[,3])*pi/180
            long_rec <- as.numeric(coord_rec[,4])*pi/180
            he_rec <- as.numeric(coord_rec[,5])
            pos_R_geo <<- matrix(c(long_rec,lat_rec,he_rec),ncol=3)
      }
      else{   # si le recepteur est statique
            pos_R_geo <<- matrix(pos_R_geo,ncol=3)
      }

      pos_R_geo_dyn <- pos_R_geo
      
      liste_f_complet <- list.files(path = path_output, pattern = (".txt"))# on va traiter le fichier present dans le dossier 2.output
            if(length(liste_f_complet)==0 | !('sortie_f.txt' %in% liste_f_complet)){
                  print('Il n\'y a aucun fichier \'sortie_f.txt\' a traiter dans le dossier 2.output',quote=FALSE); 
            }
      
      # on recupere les infos du fichier de sortie de f
      mon_data_f <-matrix(scan(paste(path_output,'sortie_f.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=14,byrow=TRUE)
      ma_date <- as.POSIXct(paste(mon_data_f[,1],'-',mon_data_f[,2],'-',mon_data_f[,3],' ',mon_data_f[,4],':',mon_data_f[,5],':',mon_data_f[,6],sep=''),format="%Y-%m-%d %H:%M:%S")
      nom_sat <- mon_data_f[,7]

      
      path_configuration <<- gsub(' ','',paste(getwd(),'/Simulateur/Configuration/')) # path of the Configuration folder
      path_input_pos_sat <<- gsub(' ','',paste(getwd(),'/Simulateur/Input/satellites_position/')) # path of the folder Input/Satellites_positions
      path_input_amfx <<- gsub(' ','',paste(getwd(),'/Simulateur/Input/AMFX/')) # Path of the folder AMFX
      path_input_mnt <<-  gsub(' ','',paste(getwd(),'/Simulateur/Input/DEM/')) # Path of the folder Input/DEM
      path_simu_programs <<- gsub(' ','',paste(getwd(),'/Simulateur/Programs/')) # path of the Configuration folder
      
      if(system_exploitation=='windows'){
             system('cmd /c copie_fichier_pos.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)      
      }
      if(system_exploitation=='linux'){
            system('cp ../0.input/*.pos Simulateur/Input/Satellites_position/')
      }
      
  
      require(tcltk)
      tclRequire("BWidget")
      
      source(paste(path_simu_programs,'simu_enreg_config.r',sep=''))
      source(paste(path_simu_programs,'pos2geo.R',sep=''))
      source(paste(path_simu_programs,'geo2pos.R',sep=''))
      source(paste(path_simu_programs,'txt2kml.R',sep=''))
      source(paste(path_simu_programs,'conversion_alti_he.R',sep=''))

            
           
      param <- matrix(seq(1,80,1),ncol=2)
      param[2,2] <- '-'    #nom_fichier_mnt
      param[3,2] <- 0      #Heure de debut simulation
      param[4,2] <- 0      # minute de debut simulation
      param[5,2] <- 23      #Heure de fin simulation
      param[6,2] <- 45     #Heure de fin simulation
      param[7,2] <- hauteur_antenne_simu  #hauteur d'antenne
      param[8,2] <- pos_R_geo[1]*180/pi    # longitude recepteur
      param[9,2] <- pos_R_geo[2]*180/pi   # latitude recepteur
      param[10,2] <- pos_R_geo[3]  # hauteur ellipsoidale recepteur
      param[11,2] <- 1e-7  #tolerance dans le calcul avec MNT
      param[12,2] <- 0     #masque elevation minimum
      param[13,2] <- 90    # masque elevation maximum
      param[14,2] <- 0      #masque azimuth minimum
      param[15,2] <- 360    # masque azimuth maximum
      param[16,2] <- 'yes'   #integrer GLONASS
      param[17,2] <- algo_simu   # algorithme
      param[18,2] <- '-'  #Free
      param[19,2] <- 1e-7  # tolerance avec ellipsoide
      param[20,2] <- 1e-7   # tolerance avec sphere
      param[21,2] <- 6378137   #semi major axis WGS84
      param[22,2] <- 298.257223563  # 1/f WGS84
      param[24,2] <- 7.292115e-5  # rotation speed of the Earth
      param[25,2] <- 398600441800000 #gravitationnal constant
      param[26,2] <- 0  # longitude du coin bas gauche du MNT
      param[27,2] <- 0  # latitude du coin bas gauche du MNT
      param[28,2] <- 0  # longitude du coin haut droit du MNT
      param[29,2] <- 0  # latitude du coin haut droit du MNT
      param[30,2] <- 1 #new sampling rate of the DEM
      param[31,2] <- 'cubic'  # mode dinterpolation MNT
      param[32,2] <- 'yes'   #integrer GPS
      param[33,2] <- 'yes'   # integrer Galileo
      param[34,2] <- 0.1902  # longueur onde
      param[35,2] <- 'Yes'  #calculer les surfaces de Fresnel
      param[36,2] <- 'No' #afficher la progression des calculs
      param[37,2] <- 'no'  #tropospheric correction
      param[38,2] <- '-'   #nom de l'AMFX calcule jusquaux sats
      param[39,2] <- '-'   #nom de l'AMFX calcule jusqua la station
      param[40,2] <- '-'   # nom de la station dans le fichier d'AMFX
      couleur_fond <<- rgb(0.88235,0.9451,0.992)

      data_simu_gps <- list()
      data_simu_glonass <- list()
      data_simu_galileo <- list()
      incr_simu <- 1
      library(lubridate,quietly=TRUE)
      mon_data_final <- c(1,2,3,4,5,6,7)
      
      for(i in 1:length(ma_date)){
            date_cur <- ma_date[i]
            nom_sat_cur <- nom_sat[i]
            date_cur <- as.POSIXct(date_cur,origin='1970-01-01')
            nom_fichier_pos_sat <- determine_sp3_correspondant(year(date_cur),month(date_cur),day(date_cur),year(date_cur),month(date_cur),day(date_cur))[1]  
            liste_sp3 <- list.files(path=path_input_pos_sat,pattern='.sp3')
            param[23,2] <- paste(day(date_cur),month(date_cur),year(date_cur),sep='/')
            param[1,2] <- paste('grg',nom_fichier_pos_sat,sep='')
            
            
            if(!(paste('grg',nom_fichier_pos_sat,sep='') %in% liste_sp3)){
                  filename <- paste(nom_fichier_pos_sat,'.Z',sep='')  
                  telecharger_sp3(filename,path_input_pos_sat)
            }
            if(system_exploitation == "windows"){
              system('cmd /c del_z_simu.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
            }
            if(system_exploitation == "linux"){
                  system('gunzip ../0.input/*.Z')
                  system('rm Simulateur/Input/Satellites_position/*.Z') 
            }
      
            if(!recepteur_statique){
                  param[1,2] <- fichier_pos
                  long_recepteur <- approx(date_rec,pos_R_geo_dyn[,1],as.double(date_cur))$y
                  lat_recepteur <- approx(date_rec,pos_R_geo_dyn[,2],as.double(date_cur))$y
                  he_recepteur <- approx(date_rec,pos_R_geo_dyn[,3],as.double(date_cur))$y
                  pos_R_geo <<- c(long_recepteur,lat_recepteur,he_recepteur)
                  pos_R_xyz <<- geo2pos(R_0_WGS84,R_p_WGS84,pos_R_geo)
                  param[8,2] <- pos_R_geo[1]*180/pi    # longitude recepteur
                  param[9,2] <- pos_R_geo[2]*180/pi   # latitude recepteur
                  param[10,2] <- pos_R_geo[3] 
            }
             
#            if(system_exploitation=='windows'){
#                  ecrit <- paste('CD Simulateur/Output/',sep='')
#                  write.table(ecrit,'simu_del_output.bat',append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)
#                  ecrit <- paste('rmdir',gsub('/','_',param[23,2]),'/S/Q')
#                  write.table(ecrit,'simu_del_output.bat',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
#                  system('cmd /c simu_del_output.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
#            }
#            if(system_exploitation=='linux'){
#                  system('rm ../Output/*')
#            }
            nom_fichier_param_temp <- gsub(' ','',paste(path_configuration,'parameters_temp.tmp'))
            simu_enreg_config(nom_fichier_param_temp,param[1,2],param[2,2],param[3,2],param[4,2],param[5,2],param[6,2],param[7,2],param[8,2],param[9,2],param[10,2],param[11,2],param[12,2],param[13,2],param[14,2],param[15,2],param[16,2],param[17,2],param[18,2],param[19,2],param[20,2],param[21,2],param[22,2],param[23,2],param[24,2],param[25,2],param[26,2],param[27,2],param[28,2],param[29,2],param[30,2],param[31,2],param[32,2],param[33,2],param[34,2],param[35,2],param[36,2],param[37,2],param[38,2],param[39,2],param[40,2])
            source(paste(path_simu_programs,"lect_sp3.r",sep=''))
            result_lect <- lect_sp3(gsub(' ','',paste(path_input_pos_sat,paste('grg',nom_fichier_pos_sat,sep=''))),param[3,2],param[4,2],param[5,2],param[6,2],param[32,2],param[16,2],param[33,2])
            nombre_satellites <- result_lect[[1]]
            mon_data <- result_lect[[2]]
            ma_date_cur <-c('*',year(date_cur),month(date_cur),day(date_cur),hour(date_cur),minute(date_cur),second(date_cur))
            

            mon_data_date <- mon_data[mon_data[,1]=='*',]
            mon_data_data <- mon_data[mon_data[,1]==paste('P',nom_sat_cur,sep=''),]
            mon_data_date <- (as.POSIXct(paste(mon_data_date[,2],'-',mon_data_date[,3],'-',mon_data_date[,4],' ',mon_data_date[,5],':',mon_data_date[,6],':',mon_data_date[,7],sep='',tz='UTC')))
            X <- approx(mon_data_date,mon_data_data[,2],(date_cur))$y
            Y <- approx(mon_data_date,mon_data_data[,3],(date_cur))$y
            Z <- approx(mon_data_date,mon_data_data[,4],(date_cur))$y
            
            mon_data <- matrix(c(ma_date_cur,paste('P',nom_sat_cur,sep=''),X,Y,Z,NA,NA,NA),ncol=7,byrow=TRUE)
            mon_data_final <- rbind(mon_data_final,mon_data)
            
      }
            mon_data <- mon_data_final[2:dim(mon_data_final)[1],]
            nb_ligne <- dim(mon_data)[1]+1  
            setwd(paste(getwd(),'/Simulateur/Output/',sep='')) 
            dir.create(gsub('/','_',param[23,2]),showWarnings=FALSE)
            setwd(paste(getwd(),'../../../',sep=''))
            path_output_ini <- path_output
            path_output <<- paste(getwd(),'/Simulateur/Output/',gsub('/','_',param[23,2]),'/',sep='') 
            setwd(paste(getwd(),'/Simulateur/Programs/',sep=''))  
            source('simu_principal.r')
            principal(nom_fichier_param_temp,mon_data,nombre_satellites,nb_ligne)
            setwd(paste(getwd(),'../../../',sep='')) 

            if(system_exploitation=='windows'){
                  ecrit <- paste('CD Simulateur/Output/',gsub('/','_',param[23,2]),'/',sep='')
                  write.table(ecrit,'simu_del_temp.bat',append=FALSE,quote=FALSE,row.names=FALSE,col.names=FALSE)
                  ecrit <- 'DEL *.tmp'
                  write.table(ecrit,'simu_del_temp.bat',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
                  ecrit <- 'COPY ..\\..\\Programs\\images\\rec.png '
                  write.table(ecrit,'simu_del_temp.bat',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
                  ecrit <- 'COPY ..\\..\\Programs\\images\\sat.png '
                  write.table(ecrit,'simu_del_temp.bat',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
                  system('cmd /c simu_del_temp.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
            }
            if(system_exploitation=='linux'){
                  system('rm ../Output/*.tmp')
                  system(paste('cp ..\\..\\Programs\\images\\rec.png ',paste(gsub('/','_',param[23,2]),'/',sep=''),sep=''))
                  system(paste('cp ..\\..\\Programs\\images\\sat.png ',paste(gsub('/','_',param[23,2]),'/',sep=''),sep=''))
            }

            data_simu_gps <- matrix(scan(paste('Simulateur/Output/',gsub('/','_',param[23,2]),'/specular_points_gps.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=20,byrow=TRUE)
            data_simu_glonass <- matrix(scan(paste('Simulateur/Output/',gsub('/','_',param[23,2]),'/specular_points_glonass.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=20,byrow=TRUE)
            data_simu_galileo <- matrix(scan(paste('Simulateur/Output/',gsub('/','_',param[23,2]),'/specular_points_galileo.txt',sep=''),what=character(0),skip=1,quiet=TRUE),ncol=20,byrow=TRUE)
            if(length(data_simu_gps)==0){
                  data_simu_gps <- matrix(rep(NaN,20),ncol=20)
            }
            if(length(data_simu_glonass)==0){
                  data_simu_glonass <- matrix(rep(NaN,20),ncol=20)
            }
            if(length(data_simu_galileo)==0){
                  data_simu_galileo <- matrix(rep(NaN,20),ncol=20)
            }
                
      
    
    
      fichier_sortie <- paste(path_output_ini,'sortie_simu.txt',sep='')
      col1 <- c('Year',format(as.numeric(year(ma_date)),digits=1))
      col2 <- c('Month',month(ma_date))
      col3 <- c('Day',day(ma_date))
      col4 <- c('Hour',hour(ma_date))
      col5 <- c('Minute',minute(ma_date))
      col6 <- c('Second',format(second(ma_date),digits=1))
      col7 <- c('PRN',nom_sat)
      col8 <- c('long(deg)',format(data_simu_gps[,1],digits=8),format(data_simu_glonass[,1],digits=8))
      col9 <- c('lat(deg)',format(data_simu_gps[,2],digits=8),format(data_simu_glonass[,2],digits=8))
      col10 <- c('he(m)',format(data_simu_gps[,3],digits=10),format(data_simu_glonass[,3],digits=10))
      col11 <- c('Semi_major_axis_Fresnel(m)',format(data_simu_gps[,9],digits=10),format(data_simu_glonass[,9],digits=8))
      col12 <- c('Semi_minor_axis_Fresnel(m)',format(data_simu_gps[,10],digits=10),format(data_simu_glonass[,10],digits=10))
      col13 <- c('Area_Fresnel(m2)',format(data_simu_gps[,11],digits=10),format(data_simu_glonass[,11],digits=10))
      col14 <- c('f',mon_data_f[,14])
      col15 <- c('Amplitude',mon_data_f[,11])
      col16 <- c('Phase',mon_data_f[,12])
      col17 <- c('Elevation(rad)',mon_data_f[,8])
      col18 <- c('Azimuth(rad)',mon_data_f[,9])
      write.matrix(matrix(c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18),ncol=18), file = fichier_sortie,sep='\t')
      
      
      closeAllConnections()
      T_fin_simu <- Sys.time()
      print('',quote=FALSE)
      print('',quote=FALSE)
      print('La simulation est finie',quote=FALSE)
      print(paste('Temps dexecution de la simulation:',difftime(T_fin_simu,T_debut_simu)),quote='FALSE')
      print('',quote=FALSE)
      print('',quote=FALSE)
      print('',quote=FALSE)

}
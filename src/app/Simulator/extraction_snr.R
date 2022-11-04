#### Programme tres tres lent mais j'arrive pas a faire mieux.


##################################################################
#
# Programme qui extrait le SNR, et calcule lelevation et lazimuth d'un rinex 3.0 + fichier sp3
#
# Input : 
#   - rinex_cur : le nom + path du fichier rinex
#   - sp3_cur, sp3_suiv, sp3_prec : le nom des fichiers d'ephemeride du jour correspondant au fichier rinex, ainsi que du jour d'avant et du jour d'apres (pour interpoler les coordonnees du satellite)
# Output :                
#   - cree un fichier extraction_[...] (dans le dossier 1.extraction) par satellite contenant toutes les infos.
#
#   juillet 2014 - Nicolas ROUSSEL
#
##################################################################

extraction <- function(rinex_cur,annee_debut,mois_debut,jour_debut,anne_fin,mois_fin,jour_fin){

        sp3 <- determine_sp3_correspondant(annee_debut,mois_debut,jour_debut,annee_fin,mois_fin,jour_fin)  # on determine le fichier sp3 correspondant a chaque fichier dobservation
        sp3_cur <- sp3[[1]] # le fichier sp3 correspondant au fichier dobservation
        sp3_prec <- sp3[[2]] # le fichier sp3 correspondant au jour d'avant le fichier d'observation (pour pouvoir faire une bonne interpolation des coordonnees satellites)
        sp3_suiv <- sp3[[3]] # le fichier sp3 correspondant au jour d'apres le fichier d'observation (pour pouvoir faire une bonne interpolation des coordonnees satellites)
        if(system_exploitation == "windows"){
              system('cmd /c del_sp3.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
              system('cmd /c del_z.bat',show.output.on.console=FALSE,ignore.stderr=TRUE,ignore.stdout=TRUE)
        }
        if(system_exploitation == "linux"){
              system('rm ../0.input/*.sp3') 
        }
	
        liste_sp3 <- list.files(path=path_input,pattern='.sp3')
      
        print('Telechargement des ephemerides',quote=FALSE)
        if(!(sp3_prec %in% liste_sp3)){
          sp3_prec<-telecharger_sp3(sp3_prec,path_input)
        }
        if(!(sp3_suiv %in% liste_sp3)){
          sp3_suiv = telecharger_sp3(sp3_suiv,path_input)
        }
        if(!(sp3_cur %in% liste_sp3)){
          sp3_cur = telecharger_sp3(sp3_cur,path_input)
        }
    
    ####
    # On lit le fichier rinex
    ####
    
    print(paste('Extraction du fichier rinex ',rinex_cur,' en cours avec l\'ephemeride ',sp3_cur,'',sep=''),quote=FALSE)  
    fichier <- file( paste(path_input,rinex_cur,sep=''), open = "r")
    lines <- readLines(fichier,warn=FALSE)
    version <- unlist(strsplit(gsub("\\s+", " ",lines[1]),' '))[2]
    
    
    ma_seq <- seq(1, length(lines),1)
    fin_header <- ma_seq[regexpr('END OF HEADER',lines)!=-1]
    header_rinex <- lines[1:fin_header]
    lines <- lines[(fin_header+1):length(lines)]
    version <- strsplit(gsub("\\s+", " ",header_rinex),' ')[[1]][2]
    if(as.numeric(version)<3){
          print('Erreur : le rinex doit etre version 3.00 au minimum',quote=FALSE);PLOUF
    }
    close(fichier)
    
    # on recupere les ent?tes d'organisation du rinex
    print(paste('Analyse du rinex ',rinex_cur,sep=''),quote=FALSE)
    a_pre <- (regexpr('OBS TYPES',header_rinex))
    ma_seq <- seq(1, length(header_rinex),1)
    a_pre <- ma_seq[a_pre!=-1]
    organisation <- header_rinex[a_pre]
    organisation_G <- 0
    organisation_R <- 0
    organisation_E <- 0
    if((length(organisation)-1)>1){
          for (u in 1:(length(organisation)-1)){
                if(substr(organisation[u],1,1)=='G'){
                        nb_info <- strsplit(gsub("\\s+", " ",organisation[u]),' ')[[1]]
                        nb_info_G <- as.numeric(nb_info[2])
                        organisation_G <- organisation[u]
                        if(!is.na(organisation[u+1])){
                              if(substr(organisation[u+1],1,1)==' '){
                                    organisation_G <- paste(organisation_G,organisation[u+1])      
                              }
                        }      
                }
                else if(substr(organisation[u],1,1)=='R'){
                      organisation_R <- organisation[u]
                      if(!is.na(organisation[u+1])){
                            if(substr(organisation[u+1],1,1)==' '){
                                  organisation_R <- paste(organisation_R,organisation[u+1])      
                            }  
                      }    
                }
                else if(substr(organisation[u],1,1)=='E'){
                      organisation_E <- organisation[u]
                      if(!is.na(organisation[u+1])){
                            if(substr(organisation[u+1],1,1)==' '){
                                  organisation_E <- paste(organisation_E,organisation[u+1])      
                            }   
                      }   
                }
          }
    }
    
 	  else{
      		if(substr(organisation,1,1)=='G'){
      			organisation_G <- organisation
      		}
      		else if(substr(organisation,1,1)=='R'){
      			organisation_R <- organisation
      		}
      		else if(substr(organisation,1,1)=='E'){
      			organisation_E <- organisation
      	 	}
   }

    if(substr(organisation[length(organisation)],1,1)=='G'){
          organisation_G <- organisation[length(organisation)]   
    }
    if(substr(organisation[length(organisation)],1,1)=='R'){
          organisation_R <- organisation[length(organisation)]   
    }
    if(substr(organisation[length(organisation)],1,1)=='E'){
          organisation_E <- organisation[length(organisation)]   
    }
    if(substr(organisation[length(organisation)],1,1)==' '){
          deb <- substr(organisation[length(organisation)-1],1,1)
          if(deb=='G'){
                organisation_G <- paste(organisation_G,organisation[length(organisation)])
          }
          else if(deb=='R'){
                organisation_R <- paste(organisation_R,organisation[length(organisation)])
          }
          else if(deb=='E'){
                organisation_E <- paste(organisation_E,organisation[length(organisation)])
          }  
    }

organisation_G <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",gsub('SYS / # / OBS TYPES','',organisation_G))),' ')[[1]]    
organisation_R <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",gsub('SYS / # / OBS TYPES','',organisation_R))),' ')[[1]]
organisation_E <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",gsub('SYS / # / OBS TYPES','',organisation_E))),' ')[[1]]


lines_date <- strsplit(gsub("\\s+", " ", gsub("^\\s+|\\s+$", "",lines[substr(lines,1,1)=='>'])),' ')
#lines_data <<- (strsplit(gsub("\\s{1,16}", " ", gsub("^\\s+|\\s+$", "",lines[(substr(lines,1,1)=='G' | substr(lines,1,1)=='R' | substr(lines,1,1)=='E' | substr(lines,1,1)=='S')])),' '))
a <- lines[(substr(lines,1,1)=='G' | substr(lines,1,1)=='R' | substr(lines,1,1)=='E' | substr(lines,1,1)=='S')]
seq_deb <- c(1,6,20,36,52,68,84,100,116,132,148,164,180,196,212,228,244,260,276,292)
seq_fin <- c(3,17,33,49,65,81,97,113,129,145,161,177,193,209,225,241,257,273,289,305)
lines_data <- matrix(rep(1,length(a)*length(seq_deb)),ncol=length(seq_deb))
for (u in 1:length(seq_deb)){
      lines_data[,u] <- substring(a,seq_deb[u],seq_fin[u])
}
lines_data <- gsub("^\\s+|\\s+$", "",lines_data,' ')
annee <- sapply(lines_date,'[[',2)
mois <- sapply(lines_date,'[[',3)
jour <- sapply(lines_date,'[[',4)
heure <- sapply(lines_date,'[[',5)
minute <- sapply(lines_date,'[[',6)
seconde <- sapply(lines_date,'[[',7)
options(digits.sec=2)
date_i <- as.double(strptime(paste(annee,'-',mois,'-',jour,' ',heure,':',minute,':',seconde,sep=''),format="%Y-%m-%d %H:%M:%OS",tz='UTC'))-17 # -17 = pour passer de GPST a UTC
nb_info <- as.numeric(sapply(lines_date,'[[',9))

ddate <- rep(date_i,nb_info)

print(paste('Analyse des ephemerides'),quote=FALSE)

if(length(nom_sat_liste_G)!=0){
      data_G <- lecture_donnees(nom_sat_liste_G,lines_data,ddate,organisation_G)
}
if(length(nom_sat_liste_R)!=0){
      data_R <- lecture_donnees(nom_sat_liste_R,lines_data,ddate,organisation_R)
}
if(length(nom_sat_liste_E)!=0){
      data_E <- lecture_donnees(nom_sat_liste_E,lines_data,ddate,organisation_E)
}

#########
# On va recuperer les coordonnees des satellites
#########
print(paste('Interpolation des coordonnees des satellites'),quote=FALSE)      
prefixe_sp3 <- 'grg'
if(length(nom_sat_liste_G)!=0){
      data_G <- calcul_coord(data_G,sp3_prec,sp3_cur,sp3_suiv)
}

if(length(nom_sat_liste_R)!=0){
      data_R <- calcul_coord(data_R,sp3_prec,sp3_cur,sp3_suiv)}
if(length(nom_sat_liste_E)!=0){
      data_E <- calcul_coord(data_E,sp3_prec,sp3_cur,sp3_suiv)
}

#########
# On calcule l'elevation et l'azimuth des satellites et on enregistre les resultats
#########


if(!recepteur_statique){
      coord_rec <- scan(nom_fichier_pos,what=character(0),quiet=TRUE)
      coord_rec <- coord_rec[(which(coord_rec=='ratio') + 1) : length(coord_rec)]
      coord_rec <- matrix(coord_rec,ncol=15,byrow=TRUE)
      date_rec <<-  as.double(as.POSIXct(paste(gsub('/','-',coord_rec[,1]),' ',coord_rec[,2],sep='')))
      lat_rec <- as.numeric(coord_rec[,3])
      long_rec <- as.numeric(coord_rec[,4])
      he_rec <- as.numeric(coord_rec[,5])
      pos_R_geo <<- matrix(c(long_rec,lat_rec,he_rec),ncol=3)
}
else{
      pos_R_geo <<- matrix(pos_R_geo,ncol=3)
}

print(paste('Calcul de lazimuth et de lelevation des satellites'),quote=FALSE)

if(length(nom_sat_liste_G)!=0){
      calcul_azimuth_elevation_et_enregistre(nom_sat_liste_G,data_G)
}
if(length(nom_sat_liste_R)!=0){
    calcul_azimuth_elevation_et_enregistre(nom_sat_liste_R,data_R)
}
if(length(nom_sat_liste_E)!=0){
      calcul_azimuth_elevation_et_enregistre(nom_sat_liste_E,data_E)
}

}
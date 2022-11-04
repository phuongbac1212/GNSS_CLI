##################################################################
#
# janvier 2016 - Nicolas ROUSSEL
#
##################################################################



directeur <<- scan(paste(path_main,'/directeur.txt',sep=''),sep='',what=character(0),quiet=TRUE,blank.lines.skip=TRUE,comment.char='#')
directeur <<- matrix(directeur,ncol=2,byrow=TRUE)

system_exploitation <<- as.character(directeur[1,2])

faire_EXTRACTION <<- switch(directeur[2,2],'oui'=TRUE,'non'=FALSE)
if(directeur[2,2]!='oui' & directeur[2,2]!='non'){print('Probleme avec le choix du programme a lancer, veuillez verifier. Valeurs possibles : \"oui\" ou \"non\"')}

faire_ANALYSE <<- switch(directeur[3,2],'oui'=TRUE,'non'=FALSE)
if(directeur[3,2]!='oui' & directeur[3,2]!='non'){print('Probleme avec le choix du programme a lancer, veuillez verifier. Valeurs possibles : \"oui\" ou \"non\"')}

faire_SIMULATEUR <<- switch(directeur[4,2],'oui'=TRUE,'non'=FALSE)
if(directeur[4,2]!='oui' & directeur[4,2]!='non'){print('Probleme avec le choix du programme a lancer, veuillez verifier. Valeurs possibles : \"oui\" ou \"non\"')}

nom_sat_liste <<- matrix(scan(paste(path_main,'liste_satellites.txt',sep=''),what=character(0),quiet=TRUE),ncol=1)[,1]
nom_sat_liste_G <<- nom_sat_liste[substr(nom_sat_liste,1,1)=='G']
nom_sat_liste_R <<- nom_sat_liste[substr(nom_sat_liste,1,1)=='R']
nom_sat_liste_E <<- nom_sat_liste[substr(nom_sat_liste,1,1)=='E']

recepteur_statique <<- switch(directeur[5,2],'oui'=TRUE,'non'=FALSE)
nom_fichier_pos <<- paste(path_input,directeur[6,2],sep='')
x_rec <- as.numeric(directeur[7,2]) # Coordonnees du recepteur (importants, car cest a partir de ces coordonnees que sont calcules lelevation et lazimuth des sats)
y_rec <-  as.numeric(directeur[8,2])
z_rec <- as.numeric(directeur[9,2])
pos_R_xyz <<- c(x_rec,y_rec,z_rec)

source('pos2geo.r')
source('geo2pos.r')
result <- pos2geo(R_0_WGS84,R_p_WGS84,pos_R_xyz)
pos_R_geo <<-  c(result[1]*pi/180,result[2]*pi/180,result[3])
url_igs <<- directeur[10,2]   # site FTP pour recuperer les ephemerides
longueur_minimum_serie <<- as.numeric(directeur[11,2])
trou_max <<- as.numeric(directeur[12,2])

path_seven_zip <<- path_Z_zip #gsub('_',' ',directeur[13,2])

methode_a_utiliser <<- directeur[13,2]
liste_code_a_utiliser <<- strsplit(directeur[14,2],',')[[1]] #*
if(length(liste_code_a_utiliser)==1){
      if(liste_code_a_utiliser=='*'){
            liste_code_a_utiliser <<- liste_codes_possibles      
      }
}
elevation_min <<- as.numeric(directeur[15,2])*pi/180
elevation_max <<- as.numeric(directeur[16,2])*pi/180
masque_azimuth <<- as.numeric(strsplit(directeur[17,2],'/')[[1]])*pi/180
h_min <<- as.numeric(directeur[18,2])
h_max <<- as.numeric(directeur[19,2])
h_point_max <<- as.numeric(directeur[20,2])
n_period <<- as.numeric(directeur[21,2])
h_fixed <<- switch(directeur[22,2],'oui'=TRUE,'non'=FALSE)
H0 <<- as.numeric(directeur[23,2])
Th<- as.numeric(directeur[24,2])*60
Fh <<- 1/Th
ecart_max <<- as.numeric(directeur[25,2])
hauteur_antenne_simu <<- as.numeric(directeur[26,2])
algo_simu <<- directeur[27,2]


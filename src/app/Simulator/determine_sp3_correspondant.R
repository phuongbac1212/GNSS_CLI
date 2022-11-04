##################################################################
#
#  Determination du fichier d'ephemeride (sp3) correspondant a un fichier d'observation RINEX
#
# Input : 
#   - fichier_rinex : le nom du fichier rinex
# Output :                
#   - sp3_correspondant : le nom du fichier d'ephemerides (sp3) du jour correspondant au fichier d'observation RINEX
#   - sp3_correspondant_prec : le nom du fichier d'ephemerides (sp3) de la veille du jour correspondant au fichier d'observation RINEX
#   - sp3_correspondant_suiv : le nom du fichier d'ephemerides (sp3) du lendemain du jour correspondant au fichier d'observation RINEX
#
#   juillet 2014 - Nicolas ROUSSEL
#
##################################################################


determine_sp3_correspondant <- function(annee_debut,mois_debut,jour_debut,annee_fin,mois_fin,jour_fin){

  date_GPS_debut <- jul_2_gps(annee_debut,mois_debut,jour_debut)
  date_GPS_fin <- jul_2_gps(annee_fin,mois_fin,jour_fin)
  #if(date_GPS_debut != date_GPS_fin){
  #  print('Il faut un fichier rinex par jour maximum. Si le fichier rinex setend sur plusieurs jours, le programme bugue');
  #  CRASH
  #}

  # if((as.numeric(date_GPS_debut) - as.numeric(date_GPS_fin))>1){
  #   print('Il faut un fichier rinex de 24h maximum. Sinon le programme bugue');
  #   CRASH
  # }

  sp3_correspondant <- paste(date_GPS_debut,'.obs',sep='')
  sp3_correspondant_prec <-  paste(jul_2_gps(annee_debut,mois_debut,jour_debut-1),'.obs',sep='')
  sp3_correspondant_suiv <-  paste(jul_2_gps(annee_debut,mois_debut,jour_debut+1),'.obs',sep='')

  return(c(sp3_correspondant,sp3_correspondant_prec,sp3_correspondant_suiv))
}

##################################################################
#
#  Determination de la longueur d'onde du signal SNR analyse
#
# Input : 
#   - maliste : vecteur contenant la liste des satellites
#   - moncode : le code qui va faire l'objet de mesure (exemple : S1C, S2W, etc...)
# Output :                
#   - longueur_onde : vecteur contenant la longueur d'onde correspondante a chaque satellite analyse et du code SNR etudie
#
#   juillet 2014 - Nicolas ROUSSEL
#
##################################################################


determine_lambda <- function(maliste,moncode){

      moncode <- substr(moncode,1,2)
      numero_de_mon_sat <- paste('_',maliste,sep='')
      n <- as.numeric(substr(numero_de_mon_sat,3,4))
      constel <- substr(numero_de_mon_sat,2,2) #'G' --> GPS, 'R' --> GLONASS

      # Pour les satellites GLONASS, on cherche le canal sur lequel ils emettent, pour savoir sur quelle frequence ils emettent
      k <- n
      k[n==1 & constel == 'R']<- 1
      k[n==2 & constel == 'R']<- -4
      k[n==3 & constel == 'R']<- 5
      k[n==4 & constel == 'R']<- 6
      k[n==5 & constel == 'R']<- 1
      k[n==6 & constel == 'R']<- -4
      k[n==7 & constel == 'R']<- 5
      k[n==8 & constel == 'R']<- 6
      k[n==9 & constel == 'R']<- -2
      k[n==10 & constel == 'R']<- -7
      k[n==11 & constel == 'R']<- 0
      k[n==12 & constel == 'R']<- -1
      k[n==13 & constel == 'R']<- -2
      k[n==14 & constel == 'R']<- -7
      k[n==15 & constel == 'R']<- 0
      k[n==16 & constel == 'R']<- -1
      k[n==17 & constel == 'R']<- 4
      k[n==18 & constel == 'R']<- -3
      k[n==19 & constel == 'R']<- 3
      k[n==20 & constel == 'R']<- 2
      k[n==21 & constel == 'R']<- 4
      k[n==22 & constel == 'R']<- -3
      k[n==23 & constel == 'R']<- 3
      k[n==24 & constel == 'R']<- 2
      
      # Du canal, on en deduit la frequence de la porteuse
      freq_L1_R <- 1602*10^6+k*0.5625*10^6
      freq_L2_R <- 1246*10^6+k*0.4375*10^6

      # De la frequence, on en deduit la longueur d'onde
      longueur_onde_L1_R <- c_lumiere/freq_L1_R
      longueur_onde_L2_R <- c_lumiere/freq_L2_R

      # on ne conserve que la longueur d'onde correspondante au code SNR en cours d'analyse
      longueur_onde <- 0
      if(moncode=='S1'){
            longueur_onde[constel=='G'] <- longueur_onde_L1_G
            longueur_onde[constel=='R'] <- longueur_onde_L1_R
      }

      if(moncode=='S2'){
            longueur_onde[constel=='G'] <- longueur_onde_L2_G
            longueur_onde[constel=='R'] <- longueur_onde_L2_R
      }

      if(moncode=='S5'){
            longueur_onde[constel=='G'] <- longueur_onde_L5_G
            longueur_onde[constel=='R'] <- NaN
      }

      return(longueur_onde)

}


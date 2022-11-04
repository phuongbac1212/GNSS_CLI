##################################################################
#
#  Determination de la frequence minimum et maximum que le SNR est sense observer en se basant sur des valeurs a priori de la hauteur min et max du recepteur, et sur une valeur de la variation maximale de cette hauteur.
#
# Input :  
#   - elevation_cur : elevation du satellite (en rad)
#   - elevation_point_cur : vitesse angulaire du satellite (en rad/s)
#   - longueur_onde : longueur d'onde (en m)
#
# Output :                 
#   - f_min : la frequence minimum
#   - f_max : la frequence maximum
#   - signe : le "signe de cette frequence" (cette frequence peut être negative)
#
#   juillet 2014 - Nicolas ROUSSEL
#
##################################################################


determine_f_min_max <- function(elevation_cur,elevation_point_cur,longueur_onde){
  
      f_min <- 99999999999999
      f_max <- 0
      f_tot <- 0
      h_point_min <- - h_point_max

      f1 <- 2/longueur_onde*(h_min+h_point_max/elevation_point_cur*tan(elevation_cur))
      f2 <- 2/longueur_onde*(h_max+h_point_max/elevation_point_cur*tan(elevation_cur))
      f3 <- 2/longueur_onde*(h_min-h_point_max/elevation_point_cur*tan(elevation_cur))
      f4 <- 2/longueur_onde*(h_max-h_point_max/elevation_point_cur*tan(elevation_cur))
      
      f_max <- max(abs(c(f1,f2,f3,f4)))
      f_min <- min(abs(c(f1,f2,f3,f4))) 
      
      if(f1 <0 & f2< 0 & f3 <0 &  f4 <0){
            signe <- -1
      }
      else{
            signe <- 1
      } 


      return(c(f_min,f_max,signe))
}
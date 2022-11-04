##################################################################
#
# Determination of the specular reflection points positions - plane approximation
# *** planar horizontal surface algorithm ***
# Launch: source("reflexion_algo_plan.R")
#
# 1st semester 2014
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
##################################################################

reflexion_algo_plan <- function(pos_E_xyz,pos_Esol_xyz,pos_R_xyz,pos_Rsol_xyz,pos_Rfictif_xyz,hauteur,distance_plan,a,b,c0,d){
  
  # pos_E_xyz : cartesian coordinates of the transmitter (satellite)
  # pos_Esol_xyz : cartesian coordinates of the projection of the satellite on the reflecting surface (nadir)
  # pos_R_xyz : cartesian coordinates of the receiver
  # pos_Rsol_xyz : cartesian coordinates of the projection of the receiver on the reflecting surface (nadir)
  # pos_Rfictif_xyz : cartesian coordinates of a fictive position of the receiver with a negative height above the ground
  # hauteur : receiver height above the reflecting surface
  # distance_plan : distance between the satellite and the plane
  # a, b, c0, d : equation of the plane
  #
  
  source('src/func/pos2geo.R')
  source('src/func/geo2pos.R')
  
  # local system coordinates
  x_Rfictif <- 0
  y_Rfictif <- -hauteur
  x_R <- 0
  y_R <- hauteur
  x_Ro <- 0
  y_Ro <- 0
  y_E <- distance_plan + hauteur
  x_E <- sqrt((vecnorm(pos_E_xyz-pos_Rsol_xyz))^2-y_E^2)
  y_S <- 0
  
  # determination of the abscissa of the specular reflection point
  x_S <- (hauteur * x_E)/(y_E + hauteur)
  
  # conversion into WGS84 cartesian coordinates
  equ_droite <- pos_Esol_xyz - pos_Rsol_xyz
  t <- x_S/vecnorm(pos_Esol_xyz-pos_Rsol_xyz)
  pos_S_xyz <- equ_droite*t+pos_Rsol_xyz
  
  # Determination of the reflection angle
  pos_E_xy <- c(x_E,y_E)
  pos_S_xy <- c(x_S,y_S)
  pos_R_xy <- c(x_R,y_R)
  
  vec_S_E <- pos_E_xy - pos_S_xy
  vec_S_R <- pos_R_xy - pos_S_xy
  theta <- asin(abs(y_E/vecnorm(c(x_S,y_S)-c(x_E,y_E))))
  
  #
  # End of the function
  #
  
  return(list(TRUE,pos_S_xyz[1],pos_S_xyz[2],pos_S_xyz[3],theta*180/pi,0,0))
  
}

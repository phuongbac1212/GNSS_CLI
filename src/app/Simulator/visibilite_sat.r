##################################################################
#
# Function which permits to see if a satellite is visible from a given station
#
# launch: source("visibilite_sat.R")
#
# 2nd semester 2012
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

sat_visible <- function(pos_E_xyz,pos_R_xyz,R_0_WGS84,R_p_WGS84,masque_vert_min,masque_vert_max,masque_hor_min,masque_hor_max,equ_plan_local){ # Function which permits to see if a satellite is visible or not 
  pos_E_xyz = matrix(pos_E_xyz, ncol = 3)
  #
  # pos_E_xyz : cartesian coordinates of the transmitter (satellite) in the WGS84 system
  # pos_R_xyz : cartesian coordinates of the receiver in the WGS84 system
  # R_0_WGS84 : semi-minor axis of the WGS84 ellipsoid
  # R_p_WGS84 : semi-major axis of the WGS84 ellipsoid
  # masque_vert_min : low vertical mask
  # masque_vert_max : up vertical mask
  # masque_hor_min : low horizontal mask 
  # masque_hor_max : up horizontal mask
  # equ_plan_local : local plane equation
  #
  
  satellite_visible <- TRUE  # Boolean which says if the satellite is visible or not

  #
  # Verification of the vertical visibility (with the satellite elevation)
  #
  
  # Determination of the satellite elevation
  distance_plan <-  dist_plan(equ_plan_local,pos_E_xyz) # We determine the distance between the satellite and its projection on the plane going through the reception point, tangent to the ellipsoid
  elevation_satellite <- asin(distance_plan/vecnorm(pos_E_xyz-pos_R_xyz))
  
  # If the satellite elevation is lower than the minimal elevation angle or upper than the maximal inclination angle with respect to the local plane.
  if(distance_plan < abs(sin(as.numeric(masque_vert_min)*pi/180)*vecnorm(pos_E_xyz-pos_R_xyz)) | distance_plan > abs(sin(as.numeric(masque_vert_max)*pi/180)*vecnorm(pos_E_xyz-pos_R_xyz)) ){
    satellite_visible <- FALSE
    return(list(satellite_visible, NA, NA))
  }

  #
  # Verification of the horizontal visibility (with the satellite azimuth)
  #
  
  # Determination of the satellite azimuth
  pos_Nord <- c(0,0,R_p_WGS84)
  pos_Nord_projete <- proj_plan(equ_plan_local,pos_Nord)
  pos_R_xyz_projete <- proj_plan(equ_plan_local,pos_R_xyz)
  pos_E_xyz_projete <- proj_plan(equ_plan_local,pos_E_xyz)
  vec_R_Nord_projete <- c(pos_Nord_projete[1]-pos_R_xyz_projete[1],pos_Nord_projete[2]-pos_R_xyz_projete[2],pos_Nord_projete[3]-pos_R_xyz_projete[3])
  vec_R_E_projete <- c(pos_E_xyz_projete[1]-pos_R_xyz_projete[1],pos_E_xyz_projete[2]-pos_R_xyz_projete[2],pos_E_xyz_projete[3]-pos_R_xyz_projete[3])
  produit_scalaire <-  vec_R_E_projete[1]*vec_R_Nord_projete[1] + vec_R_E_projete[2]*vec_R_Nord_projete[2] + vec_R_E_projete[3]*vec_R_Nord_projete[3] # xx' + yy' +zz'
  cos_azimuth <- produit_scalaire/(vecnorm(vec_R_Nord_projete)*vecnorm(vec_R_E_projete))
  produit_vectoriel <- c(vec_R_E_projete[2]*vec_R_Nord_projete[3] - vec_R_E_projete[3]*vec_R_Nord_projete[2],vec_R_E_projete[3]*vec_R_Nord_projete[1] - vec_R_E_projete[1]*vec_R_Nord_projete[3],vec_R_E_projete[1]*vec_R_Nord_projete[2] - vec_R_E_projete[2]*vec_R_Nord_projete[1])
  distance <- dist_plan(equ_plan_local,produit_vectoriel)
  signe <- sign(distance)
  azimuth <- acos(cos_azimuth)
  if(signe < 0) {azimuth <- 2*pi - azimuth}
  
  # Verificatino of the horizontal visibility
  if(as.numeric(masque_hor_min) < as.numeric(masque_hor_max)){
    if(azimuth <= as.numeric(masque_hor_min)*pi/180 | azimuth >= as.numeric(masque_hor_max)*pi/180){
      satellite_visible <- FALSE
    }
  }
  else{
    if((azimuth >= as.numeric(masque_hor_max)*pi/180) & azimuth <= as.numeric(masque_hor_min)*pi/180){
      satellite_visible <- FALSE
    }
  }   

result <- list(satellite_visible,azimuth,elevation_satellite)
return(result)

}
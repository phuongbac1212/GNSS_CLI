##################################################################
#
# This file contain function that calculate azimuth and elevation
# as well as distance from reciever and satellite coordinate
#
# Nov 2020 - Nguyen Phuong Bac
#
##################################################################

################# DEPENDANCES ##############
source("src/func/geo2pos.R")
source("src/func/pos2geo.R")

require(geosphere, quietly=TRUE)
require(pracma, quietly=TRUE)


#' @function: to calculate azimuth and elevation from satellite and reciever coordinates in ECEF
#' @param:    x, y, z: satellite coodinates
#' @param:    xr, yr, zr: reciever coodinates
#' @param:        R_p, R_0: semi=mojor and semi-minor axis of elipsoid
#'           f ....
#' @return:   cbind(azimuth, elevation)
#' @usage: source("calculate_azimuth_elevation.R")
calculate_azimuth_elevation <-
  function(x,
           y,
           z,
           xr = rec.x,
           yr = rec.y,
           zr = rec.z,
           R_0_WGS84 = 6378137,
           R_p_WGS84 = 6356752.314245179497563967,
           f = 1 / 298.257223563) {
    ########### CALCULATE AZIMUTH ##################
    pos_R_xyz <<- c(xr, yr, zr)
    result <- pos2geo(pos_R_xyz, R_0_WGS84, R_p_WGS84)
    pos_R_geo <<-
      c(result[1] * pi / 180, result[2] * pi / 180, result[3])
    pos_R_geo_cur <-
      cbind(pos_R_geo[1] * 180 / pi, pos_R_geo[2] * 180 / pi, pos_R_geo[3]) # = result
    
    pos_R_xyz_cur <- geo2pos(pos_R_geo_cur, R_0_WGS84,R_p_WGS84) 
    pos_R_0_xyz <- geo2pos((pos_R_geo_cur+c(rep(0,dim(pos_R_geo_cur)[1]),rep(0,dim(pos_R_geo_cur)[1]),rep(-5,dim(pos_R_geo_cur)[1]))), R_0_WGS84,R_p_WGS84)
    vec_rec_0 <- rbind(pos_R_0_xyz - pos_R_xyz_cur)
    X_sat <- x
    Y_sat <- y
    Z_sat <- z
    pos_sat_xyz <- cbind(X_sat,Y_sat,Z_sat)
    coord_sphere_sat <- cart2sph(pos_sat_xyz)*180/pi
    coord_sphere_rec <- cart2sph(pos_R_xyz_cur)*180/pi
    coord_sphere_rec <<- coord_sphere_rec
    coord_sphere_rec <- matrix(coord_sphere_rec,ncol=3)
    
    azimuth <- (bearing(cbind(coord_sphere_rec[,1],coord_sphere_rec[,2]),cbind(coord_sphere_sat[,1],coord_sphere_sat[,2])))*pi/180 
    azimuth[azimuth<0] <- azimuth[azimuth<0] + 2*pi
      
    vec_rec_sat <- cbind(X_sat-pos_R_xyz_cur[,1],Y_sat-pos_R_xyz_cur[,2],Z_sat-pos_R_xyz_cur[,3])
    scalaire <- vec_rec_sat[,1]*vec_rec_0[,1] + vec_rec_sat[,2]*vec_rec_0[,2] + vec_rec_sat[,3]*vec_rec_0[,3]
    
    cos_elevation_1 <- scalaire/sqrt(vec_rec_0[,3]^2+vec_rec_0[,2]^2+vec_rec_0[,1]^2)
    cos_elevation_2 <- cos_elevation_1/sqrt(vec_rec_sat[,3]^2+vec_rec_sat[,2]^2+vec_rec_sat[,1]^2)
    acos_elevation <- acos(cos_elevation_2)
    elevation <- acos_elevation - pi/2
    elevation[elevation>=90*pi/180] <- 180*pi/180-elevation[elevation>=90*pi/180]
    
    
    # on calcule la distance sat <--> rec
    
    distance <- sqrt(apply(vec_rec_sat^2,1,sum))
    
    
    return(data.frame(
      "azimuth" = azimuth,
      "elevation" = elevation,
      "distance" = distance
    ))
  }

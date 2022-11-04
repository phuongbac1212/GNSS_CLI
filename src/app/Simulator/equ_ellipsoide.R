##################################################################
#
# Calculation of the equation of the ellipsoid "parallel" to the WGS84 ellipsoid, passing through the projection to the ground of a given point
#
# Launch: source("equ_ellipsoide.R")
#
# 2nd semester 2012
#
# Author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

equ_ellipsoide <- function(R_0_WGS84,R_p_WGS84,pos_xyz,hauteur){

#
# R_0_WGS84 : semi-major axis of the WGS84 ellipsoid
# R_p_WGS84 : semi-minor axis of the WGS84 ellipsoid
# pos_xyz : cartesian coordinates of the given point
# hauteur : height of the point above the ground
#



pos_geo <- pos2geo(R_0_WGS84,R_p_WGS84,pos_xyz) # geographic coordinates of the given point


h <- pos_geo[3] - hauteur
c <- sqrt(R_0_WGS84^2*cos(pos_geo[2]*pi/180)*cos(pos_geo[2]*pi/180)+R_p_WGS84^2*sin(pos_geo[2]*pi/180)*sin(pos_geo[2]*pi/180))
R_0 <- sqrt(R_0_WGS84^2+h^2+h*c+R_0_WGS84^2*h/c)
R_p <- sqrt(R_p_WGS84^2+h^2+h*c+R_p_WGS84^2*h/c)


#
# End of the function
# 

return (c(R_0,R_p))
}
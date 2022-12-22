##################################################################
#
# Determination of the equation of the tangent plane to the ellipsoid passing through a given point, at that given point
#
# Launch: source("equ_plan.R")
#
# 2nd semester 2012
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

equ_plan <- function(R_0_init=6378137,R_p_init=6356752.3142,pos_xyz){

#
# R_0_WGS84 : semi-major axis of the WGS84 ellipsoid
# R_p_WGS84 : semi-minor axis of the WGS84 ellipsoid
# pos_xyz : cartesian coordinates of the given point
#

# source("equ_ellipsoide.R")
# source("pos2geo.R")
# source("geo2pos.R")

#
# Calculation of the local ellipsoid : the ellipsoid parallel to the WGS84, passing through the given point
#

ellipsoide_recepteur <- equ_ellipsoide(R_0_init,R_p_init,pos_xyz,0) # equation of the local ellipsoid
R_0 <- ellipsoide_recepteur[1] # semi-major axis of the local ellipsoid
R_p <- ellipsoide_recepteur[2] # semi-minor axis of the local ellipsoid
e <-  sqrt((R_0^2-R_p^2)/R_0^2) # excentricity of the local ellipsoid

#
# Calculation of the vector normal to the plane tangent to the local ellipsoid, at the given point
# Let that vector be the vector separating the given point with the intersection of the major-normal N of the given point with the abscissa axis
#

pos_geo_ellipsoide_local <- pos2geo(R_0,R_p,pos_xyz) # geographic coordinates of the given point, in the local ellipsoid ! (and not in the WGS84 ellipsoid !)
N <- R_0/sqrt(1-e^2*sin(pos_geo_ellipsoide_local[2])*sin(pos_geo_ellipsoide_local[2])) # major-normal N of the ellipse at the given point
n <- e^2*N
pos_proj_sur_plan_equateur_xyz <- geo2pos(R_0,R_p,c(pos_geo_ellipsoide_local[1],pos_geo_ellipsoide_local[2],-(N-n))) # coordinates of the intersection between the major-normal N and the abscissa axis
coord_vec_normal_xyz <- pos_xyz - pos_proj_sur_plan_equateur_xyz # cartesian coordinates of the normal vector

#
# Equation of a plane : ax + by + cz + d = 0
#

a <- coord_vec_normal_xyz[1]
b <- coord_vec_normal_xyz[2]
c <- coord_vec_normal_xyz[3]

#
# Calculation of d, knowing that the plane passes through the given point
#

d <- -(a *pos_xyz[1] + b * pos_xyz[2] + c * pos_xyz[3]) 

#
# End of the function
#

return(c(a,b,c,d))
}
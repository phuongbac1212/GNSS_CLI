##################################################################
#
# Calculation of the coordinates of the intersection points of two straight lines tangent to a given ellipsoid and passing through a given point. 
#
# Launch: source("coord_pts_tangente.R")
#
# 2nd semester 2012
#
# Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################


coord_pts_tangente <- function(R_0,R_p,pos_R_xyz,pos_Rsol_xyz,pos_Rsol_geo,pos_R_geo){

# R_0 : semi-major axis of the ellipsoid
# R_p : semi-minor axis of the ellipsoid
# pos_R_xyz : cartesian coordinates of the receiver
# pos_Rsol_xyz : cartesian coordinates of the receiver
# pos_Rsol_geo : cartesian coordinates of the receiver with a zero altitude above the ground
# pos_R_geo : geographic coordinates of the receiver

source('tangente_ellipse.R')

x <- pos_R_xyz[1]
y <- pos_R_xyz[2]
z <- pos_R_xyz[3]

y0 <- z # ordinates of the receiver in a local system concentric with the WGS84 cartesian system and with the axis y // to Z_WGS84 and with the x axis in the meridian plane

dist_plani <- sqrt(x^2+y^2) # planimetric distance, in the (X_WGS84 - Y_WGS84) plane, between the system origin and the projection of the receiver point in the (X_WGS84 - Y_WGS84) plane
x0 <- dist_plani # abscissa of the receiver in the local system

#
# We have now the coordinates of the receiver in a local 2D system defined in the meridian plane
# We now calculate the 2 intersection points of the meridian ellips with the two straight lines tangent with this ellipse, and passing through the receiver point.
#

resultat <- tangente_ellipse(R_0,R_p,x0,y0) # we get (a1,b1,a2,b2), the parameters of the two straight lines tangent to the meridian ellipse passing throughthe receiver points, having for equation: x + a*y + b = 0
a1 <- resultat[1]
a2 <- resultat[3]
b1 <- resultat[2]
b2 <- resultat[4]

#
# We now have the equation of the two straight lines, and also the ellipse equation. We just have now to solve a system to determine the coordinates of the two points of tangence
# We have the following system: (*) x^2 / R_0^2 + y^2 / R_p^2 = 1 et (**) x + V * Y + W = 0
# That is the same to solve the following equation: y^2 * (V^2/R_0^2 + 1/b^2) + y * (2*V*W/R_0^2) + W^2/R_0^2 - 1 = 0
#

y1 <- -(2*a1*b1/R_0^2)/(2*(a1^2/R_0^2+1/R_p^2)) # ordinate of the first intersection point
y2 <- -(2*a2*b2/R_0^2)/(2*(a2^2/R_0^2+1/R_p^2)) # ordinate of the second intersection point
x1 <- -a1*y1-b1 # We deduce the abscissa of the first intersection point
x2 <- -a2*y2-b2 # We deduce the abscissa of the second intersection point

#
# End of the function
#

return(c(x1,y1,x2,y2))

}
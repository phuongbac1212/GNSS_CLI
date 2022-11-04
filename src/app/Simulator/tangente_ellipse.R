##################################################################
#
# Determination of the equation of the 2 straight lines tangent to a given ellipse and going through a given point.
# Equation looks like : x + a * y + b = 0
#
# Launch: source("reflexion.R")
#
# 2nd semester
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

tangente_ellipse <- function(R_0,R_p,x0,y0) {

#
# R_0 : semi-minor axis of the given ellipsoid
# R_p : semi-major axis of the given ellipsoid
# x0 : abscissa of the considered point in the meridian plane
# y0 : ordonate of the considered point in the meridian plane
#

#
# An ellipse equation is : x^2/R_0^2 + y^2/R_p^2 = 1 and a straight line equation is : u*x + v*y + w = 0 
# We can easily deduce that "the straight line is tangent to the ellipse if and only if R_0^2*u^2 + R_p^2*v^2 = w^2"
# We thus have the equation system : (*) R_0^2 + R_p^2 * V^2 = W^2 and (**) x0 + V*y0 + W = 0, with V = v/u and W = w/u and (x0,y0)  the coordinates of the reception point where goes through the 2 straight lines 
# We thus get the following equation, by isolating V : V^2(R_p^2-y0^2) + V(-2*x0*y0) + (R_0^2-x0^2) = 0
#

delta <- (2*x0*y0)^2 -4*(R_p^2-y0^2)*(R_0^2-x0^2)
a1 <- (2*x0*y0+sqrt(delta))/(2*(R_p^2-y0^2)) # We get a first value of V
a2 <- (2*x0*y0-sqrt(delta))/(2*(R_p^2-y0^2)) # We get a second value of V
b1 <- -x0-y0*a1 # We deduce a first value of W
b2 <- -x0-y0*a2 # We deduce a second value of W

#
# End of the function
#

return(c(a1,b1,a2,b2))

}
##################################################################
#
# Determination of the algebric distance between a point and a plane
#
# Launch: source("dist_plan.R")
#
# 2nd semester 2012
#
# Author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

dist_plan <- function(equ_plan,coord_point_xyz){
  
  #
  # Calculation of the algebric distance
  #
  
  distance <- (equ_plan[1] * coord_point_xyz[1] + equ_plan[2] * coord_point_xyz[2] + equ_plan[3] * coord_point_xyz[3] + equ_plan[4])/sqrt(equ_plan[1]^2 + equ_plan[2]^2 + equ_plan[3]^2)
  
  #
  # End of the function
  #
  
  return(distance)
  
}
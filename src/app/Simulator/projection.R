##################################################################
#
# Point projection on a given plane or ellipsoid
#
# launch: source("projection.R")
#
# 2nd semester 2012
#
# Author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
#
##################################################################

#source('mnt_fonctions.R')
proj_plan <- function(equ_plan,pos_xyz){ # point projection (cartesian coordinates pos_xyz) on the plane (equ_plan)
      # point coordinates
      xp <- pos_xyz[1]
      yp <- pos_xyz[2]
      zp <- pos_xyz[3]
      
      # plane coordinates
      a <- equ_plan[1]
      b <- equ_plan[2]
      c1 <- equ_plan[3]
      d <- equ_plan[4]
      
      # the parametric equation of the straight line orthogonal to the plane, going through P, looks like  x = a * t + xp ; y = b * t + yp ; z = z * t + zp so, by putting this parametric equation in the plane equation, we can get t:
      
      t <- -(d + xp*a + yp*b + zp*c1) / (a^2+b^2+c1^2)
      
      x <- t*a + xp
      y <- t*b + yp
      z <- t*c1 + zp
      
      return(c(x,y,z))
}

proj_surf <- function(R_0,R_p,pos_xyz){ # Projection of the point (pos_xyz) on the ellipsoid (R_0,R_p)
      
      #
      # R_0 : semi-major axis of the considered ellipsoid
      # R_p : semi-minor axis of the considered ellipsoid
      # pos_xyz : cartesian coordinates of the considered point in the considered ellipsoid
      #
      
      #source('pos2geo.R')
      #source('geo2pos.R')
      
      pos_geo <- pos2geo(R_0,R_p,pos_xyz) # geographic coordinates of the considered point in the considered ellipsoid
      pos_proj_xyz <- geo2pos(R_0,R_p,c(pos_geo[1],pos_geo[2],0)) # cartesian coordinates of the considered point projected on the considered ellipsoid (ellipsoid height set to 0)
      
      return(pos_proj_xyz)
      
      #
      # End of the function
      #

}

proj_surf_mnt <- function(R_0,R_p,pos_xyz,mnt,M_mnt,L_mnt){ # Projection of the point (pos_xyz) on a DEM

      pos_geo <- pos2geo(R_0,R_p,pos_xyz) # geographic coordinates of the point in the considered ellipsoid
      he <- interpolation(mnt,M_mnt,L_mnt,pos_geo,5) # we get the ellipsoid height of the point by interpolation with the 5 closest neighbours
      pos_proj_xyz <- geo2pos(R_0,R_p,c(pos_geo[1],pos_geo[2],he))
      
      return(pos_proj_xyz)
            
      #
      # End of the function
      #
}
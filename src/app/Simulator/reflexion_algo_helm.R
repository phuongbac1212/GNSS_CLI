##################################################################
#
# Determination of the specular reflection points positions - sphere approximation
# *** sphere algorithm ***  --> based on Helm work: Helm A.: Ground vased GPS altimetry with the L1 openGPS receiver using carrier phase-delay observations of reflected GPS signals. Thesis (PhD). GFZ. 2008
#
# Launch: source("reflexion_algo_helm.R")
#
# 2nd semester 2012
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr

reflexion_algo_helm <- function(R_T,pos_centre_xyz,pos_E_xyz,pos_R_xyz,elevation,hauteur,tolerance_helm){
      #
      # R_T : radius of the osculating sphere (gaussian radius of curvature)
      # pos_centre_xyz : cartesian coordinates of the center of the osculating sphere
      # pos_E_xyz : cartesian coordinates of the transmitter (satellite)
      # pos_R_xyz : cartesian coordinates of the receiver
      # elevation : satellite elevation (rad)
      # hauteur : receiver height (w.r.t the ground)
      # tolerance_helm : tolerance used for the calculation convergence
      #
      
      # The goal is to find gamma, the 'latitude' of the reflection point. (Well, it's not exactly the latitude. It is the angle between the equator plane and the reflection point in the plane formed by the satellite, the receiver and the center of the Earth.
      # We approximate the Earth as a sphere
      
      # Determination of the local Earth radius ==> we approximate the Earth as a sphere but we calculate a local radius to be more exact.
      
      #
      # We are in the plane P defined by the reflection point, the satellite, the receiver and the center of the Earth ==> We define a local system with the ordonate axis being the vector defined by the center of the Earth and the receiver
      # and the abscissa axis being the vector defined by the center of the Earth and a point located in the plane P such as the abscissa axis being perpendicular to the ordinate axis.
      #
      
      # Transmitter (satellite) coordinates in the local system
      pos_E_xyz = matrix(pos_E_xyz, ncol = 3)
      R_e <- vecnorm(pos_E_xyz-pos_centre_xyz)
      
      x_e <- R_e*cos(elevation)*sqrt(1-(R_T+hauteur)^2/R_e^2*cos(elevation)*cos(elevation))-(R_T + hauteur)*sin(elevation)*cos(elevation)
      y_e <- R_e*sin(elevation)*sqrt(1-(R_T+hauteur)^2/R_e^2*cos(elevation)*cos(elevation))+(R_T + hauteur)*cos(elevation)*cos(elevation)
      
      # Receiver coordinates in the local system
      x_r <- 0
      y_r <- vecnorm(pos_R_xyz-pos_centre_xyz)
      
      # angle transmitter - center of the Earth - receiver 
      beta_r_e <- acos((x_e*x_r + y_e*y_r)/(sqrt(x_e^2+y_e^2)*sqrt(x_r^2+y_r^2))) # we get the angle by dividing the scalar product by the multiplication of the norms of the two vectors
      
      # We can show that gamma checks a spherical mirror equation : c4 * t^4 + c3 * t^3 + c2 * t^2 + c1 * t + c0 = 0, with t = tan(gamma/2)
      # We can get t, and then gamma by resolving this equation.
      
      #
      # Resolution of the spherical mirror equation to get t
      #
      
      # gamma zero --> approximation of an initial gamma to fasten the convergence of the calculation 
      gamma_n <- pi/2 - 1/3*beta_r_e
      
      iterations <- 1
      correction_n <- 10000 # will by after defined by  |gamma_n_plus_1 - gamma_n|
      converge <- TRUE # says if calculation have converged or not
      
      while (correction_n > as.numeric(tolerance_helm)){ # while we did not converge sufficiently
      
        # coefficient of the spherical mirror equation c4 * t^4 + c3 * t^3 + c2 * t^2 + c1 * t + c0 = 0, with t = tan(gamma/2)
        c0 <- (x_e*y_r + y_e*x_r) - R_T*(y_e + y_r)
        c1 <- -4*(x_e*x_r - y_e*y_r) + 2*R_T*(x_e + x_r)
        c2 <- -6*(x_e*y_r + y_e*x_r)
        c3 <- 4*(x_e*x_r - y_e*y_r) + 2*R_T*(x_e + x_r)
        c4 <- (x_e*y_r + y_e*x_r) + R_T*(y_e + y_r)
        
        t_n <- tan(gamma_n/2)
        f_n <- c4*t_n^4 + c3*t_n^3 + c2*t_n^2 + c1*t_n + c0
        f_n_prime <- 4*c4*t_n^3 + 3*c3*t_n^2 + 2*c2*t_n + c1
        f_n_seconde <- 12*c4*t_n^2 + 6*c3*t_n +2*c2
      
        if (f_n_prime == 0){
          converge <- FALSE
          break
        }
      
        t_n_plus_1 <- t_n - f_n/f_n_prime
        gamma_n_plus_1 <- 2*atan(t_n_plus_1)
        correction_n <- abs(gamma_n_plus_1 - gamma_n)
        
        gamma_n <- gamma_n_plus_1
        
       	if(iterations > 50000){
       	  print('error convergence sphere approximation')
       	  converge <- FALSE
      		break
      	}	
      	iterations <- iterations + 1
      }
      
      if (!converge) {gamma_n <- 0}
      
      #
      # Coordinates of the reflection point 
      #
      
      # Coordinates of the reflection point in the local system
      x_s <- R_T*cos(gamma_n)
      y_s <- R_T*sin(gamma_n)
      arc <- (pi/2-gamma_n)*R_T  # arc length
      
      
      # Coordinates of the reflection point in the WGS84 system
      # 
      # We get :
      # X1 = A1*x1 + B1*y1 -Tx
      # Y1 = A2*x1 + B2*y1 +
      # Z1 = A3*x1 + B3*y1
      # X2 = A1*x2 + B1*y2
      # Y2 = A2*x2 + B2*y2
      # Z2 = A3*x2 + B3*y2
      # 
      # By getting A1, A2, A3, B1, B2, B3, we will be able to go from the local system to the WGS84 system
      #
      
      X1 <- pos_E_xyz[1] - pos_centre_xyz[1]
      Y1 <- pos_E_xyz[2] - pos_centre_xyz[2]
      Z1 <- pos_E_xyz[3] - pos_centre_xyz[3]
      x1 <- x_e
      y1 <- y_e
      X2 <- pos_R_xyz[1] - pos_centre_xyz[1]
      Y2 <- pos_R_xyz[2] - pos_centre_xyz[2]
      Z2 <- pos_R_xyz[3] - pos_centre_xyz[3]
      x2 <- x_r
      y2 <- y_r
      
      mat_A <- matrix(c(x1,0,0,y1,0,0,0,x1,0,0,y1,0,0,0,x1,0,0,y1,x2,0,0,y2,0,0,0,x2,0,0,y2,0,0,0,x2,0,0,y2),ncol=6,byrow=TRUE)
      mat_B <- c(X1,Y1,Z1,X2,Y2,Z2)
      mat_P <- diag(6)
      mat_X <- (solve(t(mat_A)%*%mat_P%*%mat_A)) %*%(t(mat_A)%*%mat_P%*%mat_B)
      pos_S_xyz <- c(mat_X[1]*x_s+mat_X[4]*y_s+pos_centre_xyz[1],mat_X[2]*x_s+mat_X[5]*y_s+pos_centre_xyz[2],mat_X[3]*x_s+mat_X[6]*y_s+pos_centre_xyz[3])
      
      #
      # Determination of the reflection angle
      #
      
      pos_E_xy <- c(x_e,y_e)
      pos_S_xy <- c(x_s,y_s)
      pos_R_xy <- c(x_r,y_r)
      
      vec_S_E <- pos_E_xy - pos_S_xy
      vec_S_R <- pos_R_xy - pos_S_xy
      theta <- asin(abs((y_e-y_r-hauteur)/vecnorm(c(x_s,y_s)-c(x_e,y_e))))# acos((vec_S_E[1]*vec_S_R[1]+vec_S_R[2]*vec_S_R[2])/(vecnorm(vec_S_E)*vecnorm(vec_S_R)))/2
      #
      # End of the function
      #
      
      return(list(converge,pos_S_xyz[1],pos_S_xyz[2],pos_S_xyz[3],theta*180/pi,0,R_T))

}

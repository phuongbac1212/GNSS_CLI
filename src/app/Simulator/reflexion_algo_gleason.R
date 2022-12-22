##################################################################
#
# Determination of the reflection points positions - ellipsoid approximation
# *** ellipsoid algorithm ***  --> based on Gleason work: Gleason S.: Remote Sensing of Ocean, Ice and Land Surfaces Using Bistatically Scattered GNSS Signals From Low Earth Orbit. Thesis (PhD.), University of Surrey, 2006.
#
# Launch: source("reflexion_algo_gleason.R")
#
# 2nd semester 2012
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################


signe2 <- function(E2S_norm, R2S_norm, S_norm) {
  a <- crossprod(S_norm, E2S_norm) / vecnorm(E2S_norm) / vecnorm(S_norm)
  b <-
    crossprod(S_norm, R2S_norm) / vecnorm(R2S_norm) / vecnorm(S_norm)
  if (a >= 1) {
    a <- 1
  }
  if (b >= 1) {
    b <- 1
  }
  if (a <= -1) {
    a <- -1
  }
  if (b <= -1) {
    b <- -1
  }
  theta1 <- 180 / pi * (acos(a)) # incidence angle
  theta2 <- 180 / pi * (acos(b)) # reflection angle
  if (theta1 > 90) {
    theta1 <- abs(180 - theta1)
  }
  if (theta2 > 90) {
    theta2 <- abs(180 - theta2)
  }
  
  monsigne <-   sign(theta1 - theta2)
  
  return(monsigne)
}

calcul_angles <- function(vec_1, vec_2, vec_3) {
  a <-
    crossprod(vec_1, vec_2) / vecnorm(vec_1) / vecnorm(vec_2) # incidence angle
  b <-
    crossprod(vec_3, vec_2) / vecnorm(vec_3) / vecnorm(vec_2) # reflection angle
  if (a >= 1) {
    a <- 1
  }
  if (b >= 1) {
    b <- 1
  }
  if (a <= -1) {
    a <- -1
  }
  if (b <= -1) {
    b <- -1
  }
  
  theta1 <- 180 / pi * (acos(a)) # incidence angle
  theta2 <- 180 / pi * (acos(b)) # reflection angle
  if (theta1 > 90) {
    theta1 <- abs(180 - theta1)
  }
  if (theta2 > 90) {
    theta2 <- abs(180 - theta2)
  }
  
  return(c(theta1, theta2))
}




reflexion_algo_gleason <-
  function(R_0,
           R_p,
           pos_E_xyz,
           pos_R_xyz,
           tolerance_gleason) {
    # pos_E_xyz = matrix(pos_E_xyz, ncol = 3)
    pos_R_xyz = c(pos_R_xyz)
    #
    # R_0 : semi-minor axis of the considered ellipsoid
    # R_p : semi-major axis of the considered ellipsoid
    # pos_E_xyz: cartesian coordinates of the transmitter (satellite)
    # pos_R_xyz : cartesian coordinates of the receiver
    # tol?rance_gleason : convergence tolerance of the reflection calculation
    #
    iterations <- 1
    converge <- FALSE
    
    #
    # a priori position of the specular reflection point
    #
    
    source('intersection_droite_ellipsoide.R')
    
    # we proceed by dichotomy
    ma_liste <- list()
    ma_liste[[1]] <- c(proj_surf(R_0, R_p, pos_R_xyz))
    ma_liste[[2]] <- c(proj_surf(R_0, R_p, pos_E_xyz))
    pos_precedent <- 0
    pos_S2 <- 12
    continuer <- TRUE
    while (continuer) {
      pos_S1 <- ma_liste[[1]]
      pos_S3 <- ma_liste[[2]]
      pos_precedent <- pos_S2
      pos_S2_non_proj <- (pos_S3 + pos_S1) / 2
      
      S_norm <- calcul_normale(pos_S1, pos_S3)
      pos_S2 <-
        intersection_droite_ellipsoide(
          R_0,
          R_p,
          S_norm[1],
          S_norm[2],
          S_norm[3],
          pos_S2_non_proj[1],
          pos_S2_non_proj[2],
          pos_S2_non_proj[3]
        )
      
      a <- pos2geo(R_0, R_p, pos_S2)
      S_norm <- pos_S2 - geo2pos(R_0, R_p, c(a[1], a[2], a[3] - 5))
      S_norm <- c(S_norm / vecnorm(S_norm))
      
      # Normalized differences
      E2S_norm1 <-
        (pos_E_xyz - pos_S1) / vecnorm(pos_E_xyz - pos_S1)
      R2S_norm1 <-
        (pos_R_xyz - pos_S1) / vecnorm(pos_R_xyz - pos_S1)
      E2S_norm2 <-
        (pos_E_xyz - pos_S2) / vecnorm(pos_E_xyz - pos_S2)
      R2S_norm2 <-
        (pos_R_xyz - pos_S2) / vecnorm(pos_R_xyz - pos_S2)
      E2S_norm3 <-
        (pos_E_xyz - pos_S3) / vecnorm(pos_E_xyz - pos_S3)
      R2S_norm3 <-
        (pos_R_xyz - pos_S3) / vecnorm(pos_R_xyz - pos_S3)
      
      signe_1 <- signe2(E2S_norm1, R2S_norm1, S_norm)
      signe_2 <- signe2(E2S_norm2, R2S_norm2, S_norm)
      signe_3 <- signe2(E2S_norm3, R2S_norm3, S_norm)
      
      if ((signe_1 == signe_3)) {
        #no reflection point is posible at all
        pos_precedent <- pos_S2
      }
      else if (signe_1 != signe_3) {
        if (signe_1 != signe_2) {
          # we keep the first half
          ma_liste[[1]] <- pos_S1
          ma_liste[[2]] <- pos_S2
        }
        else{
          # we keep the second half
          ma_liste[[1]] <- pos_S2
          ma_liste[[2]] <- pos_S3
        }
      }
      
      result <- calcul_angles(E2S_norm2, S_norm, R2S_norm2)
      theta1 <- result[1]
      theta2 <- result[2]
      theta_diff <- abs(theta1 - theta2)
      
      if (theta_diff < tolerance_gleason) {
        converge <-
          TRUE
        continuer <-
          FALSE
      }   # we stop convergence when the difference between incident and reflecting angle are below the tolerance.
      if (identical(pos_S2, pos_precedent)) {
        # if the dichotomious method cannot converge (possible due to the normals who are different at each of the 3 points), we proceed by iteration until reaching the position of the specular reflection point.
        # Declaration of the variables
        theta_diff <- 1200
        K <- 10
        pos_precedent <- 0
        iter <- 1
        d <- 0
        dXYZ <- 1
        theta_diff_prec <- 0
        signe <- 1
        sens_prec <- 0
        
        while (theta_diff > tolerance_gleason) {
          # while no convergence
          
          # Normalized differences
          E2S_norm <-
            (pos_E_xyz - pos_S2) / vecnorm(pos_E_xyz - pos_S2)
          R2S_norm <-
            (pos_R_xyz - pos_S2) / vecnorm(pos_R_xyz - pos_S2)
          if (iter > 1) {
            sens_prec <- sens
          }
          a <- pos2geo(R_0, R_p, pos_S2)
          S_norm <-
            pos_S2 - c(geo2pos(R_0, R_p, c(a[1], a[2], a[3] - 5)))
          S_norm <- S_norm / vecnorm(S_norm)
          
          result <-
            calcul_angles(vec_1 = E2S_norm2, vec_2 = S_norm, vec_3 = R2S_norm2)
          theta1 <- result[1]
          theta2 <- result[2]
          theta_diff <- abs(theta1 - theta2)
          theta_diff <- abs(theta1 - theta2)
          
          if (theta1 < theta2) {
            dXYZ <- R2S_norm / K
            sens <- 1
          }
          else{
            dXYZ <- E2S_norm / K
            sens <- -1
          }
          if (iter > 1) {
            if (sens_prec != sens) {
              K <- K * 10
              dXYZ <- dXYZ / K
            }
          }
          
          # Projection on the surface
          pos_S2_non_projete <-  pos_S2 + dXYZ
          
          pos_S2 <-
            c(proj_surf(R_0, R_p, pos_S2_non_projete))
          a <- pos2geo(R_0, R_p, pos_S2)
          S_norm <-
            pos_S2 - c(geo2pos(R_0, R_p, c(a[1], a[2], a[3] - 5)))
          S_norm <- S_norm / vecnorm(S_norm)
          result <-
            calcul_angles(E2S_norm2, S_norm, R2S_norm2)
          theta1 <- result[1]
          theta2 <- result[2]
          theta_diff <- abs(theta1 - theta2)
          
          if (iter > 10000) {
            print('there is an error in the ellipsoid algorithm... no convergence possible.[[1]]')
            converge <- FALSE
            break
          }
          iter <- iter + 1
        }
        if(converge == F) break
        converge <- TRUE
        continuer <- FALSE
      }
      if (iterations > 10000) {
        print('there is an error in the ellipsoid algorithm... no convergence possible.[[2]]')
        converge <- FALSE
        break
      }
      
      iterations <- iterations + 1
    }
    pos_S_xyz <- pos_S2
    
    #
    # End of the function
    #
    if (!converge) {
      print('ellipsoide pas converge')
    }

    return(list(
      converge,
      pos_S_xyz[1],
      pos_S_xyz[2],
      pos_S_xyz[3],
      mean(theta1, theta2),
      theta_diff
    ))
    
    
  }

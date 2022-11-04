##################################################################
#
# Set of functions used during the treatment of the DEM in 1D
#
# Launch: source("mnt_fonctions_1D.R")
#
# 2nd semester 2012
#
# Author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################



  projection_plan <- function(point,a,b,c,d){
      k <- (-d-a*point[1]-b*point[2]-c*point[3])/(a*a+b*b+c*c) 
      result <- c(point[1]+a*k,point[2]+b*k, point[3]+c*k)
      return(result)  
  }
  
  
  mnt_visible <- function(x,y,z,pt1,pt2,R_0_sol,R_p_sol){ # Function to know if a point 1 is visible from a point 2 or if it is hidden by the DEM
          source('calcul_normale.r')
          normale <- calcul_normale(pt1,pt2)
          a <- normale[1]
          b <- normale[2]
          c <- normale[3]
          d <- -a*pt1[1]-b*pt1[2]-c*pt1[3]
      
          
          vec_d1 <- (pt2-pt1)/vecnorm(pt2-pt1) 
          x0_d1 <- pt1

            visible <- TRUE      
            longueur <- length(x)-1
            truc <- 1
            source('equ_plan.r')

            plan_local <- equ_plan(R_0_sol,R_p_sol,pt2)

        while(visible & truc<(longueur-1)){                 
        for(num in 1:(longueur-1)){ # For each point of the DEM
          # We calculate the straight line D passing through the considered point of the DEM with the following point
              truc <- num
              x0_d2 <- c(x[num],y[num],z[num])
              x0_d2_2 <-   c(x[num+1],y[num+1],z[num+1])
          x0_d2 <- projection_plan(x0_d2,a,b,c,d)
           x0_d2_2 <- projection_plan(x0_d2_2,a,b,c,d)
          
          vec_d2 <- (x0_d2_2-x0_d2)/vecnorm(x0_d2_2-x0_d2)
          x0_d2_proj <- projection_plan(x0_d2,plan_local[1],plan_local[2],plan_local[3],plan_local[4])
          x0_d2_2_proj <- projection_plan(x0_d2_2,plan_local[1],plan_local[2],plan_local[3],plan_local[4])
           pt1_proj <- projection_plan(pt1,plan_local[1],plan_local[2],plan_local[3],plan_local[4])
           pt2_proj <- projection_plan(pt2,plan_local[1],plan_local[2],plan_local[3],plan_local[4])
         dist_ref_1 <- vecnorm(pt1_proj-x0_d2_proj)
         dist_ref_2 <- vecnorm(pt1_proj-x0_d2_2_proj)
         dist_ref_3 <- vecnorm(pt1_proj-pt2_proj)
         if(dist_ref_1<dist_ref_3 &dist_ref_2<dist_ref_3 ){
          t_prim <- (x0_d2[2]-x0_d1[2]+vec_d1[2]*x0_d1[1]/vec_d1[1]-x0_d2[1]*vec_d1[2]/vec_d1[1])/(vec_d1[2]*vec_d2[1]/vec_d1[1]-vec_d2[2])
          t <- (x0_d2[1]+vec_d2[1]*t_prim-x0_d1[1])/(vec_d1[1])
      
            pt<- c(x0_d1[1]+vec_d1[1]*t,x0_d1[2]+vec_d1[2]*t,x0_d1[3]+vec_d1[3]*t)
            dist1 <- vecnorm(pt-x0_d2)
            dist2 <- vecnorm(pt-x0_d2_2)
            dist_ref <- vecnorm(x0_d2_2-x0_d2)
   
            if(max(dist1,dist2) <= dist_ref ){
      
            visible <- FALSE
            truc <- longueur * 123132131
            break
            }
            else{
              visible <- TRUE
            }
          
      
          } 
        }
      }
      
        return(visible)
 }








distance_entre_points <- function(point_1,point_2){   # Function which gives the distance between two points
  distance <- sqrt((point_1[1]-point_2[1])^2+(point_1[2]-point_2[2])^2)
  if(distance==0) {distance <- 0.0000000000000000001}
  return(distance)
}
















intersection_droite_cercle  <- function(a,b,xi,yi, R){ # Function which calculate the two intersection point of a straight line with a circle
  delta <- (-2*xi+2*a*b-2*a*yi)^2 - 4*(1+a^2)*(xi^2+yi^2-R^2-2*b*yi+b^2)
  x1 <- (- (-2*xi + 2*a*b-2*a*yi) + sqrt(delta)) / (2*(1+a^2))
  x2 <- (- (-2*xi + 2*a*b-2*a*yi) - sqrt(delta)) / (2*(1+a^2))
  y1 <- a*x1+b
  y2 <- a*x2+b
  return(c(x1,y1,x2,y2))
}








equation_droite <- function(point_1,point_2){ # Function which calculate the equation of a straight line going through two known points
  x1 <- point_1[1]
  y1 <- point_1[2]
  x2 <- point_2[1]
  y2 <- point_2[2]
  a <- (y2 - y1)/(x2 - x1)
  b <- y1 - a * x1
  return(c(a,b))
}








interp_mnt <- function(long,lat,long_mnt,lat_mnt,he_mnt,algo){ # Function to interpolate the altitude of a point from its neighbours
library('akima')
    
 if(algo == 'bilinear'){
      alti <- interp(x=long_mnt,y=lat_mnt,z=he_mnt,xo=long,yo=lat,linear=TRUE)$z
  }
  
  else if(algo=='cubic'){
         vec_long <- 0
         vec_lat <- 0
         for(u in 1:length(long_mnt)){
            if(u!=1 &long_mnt[u]==long_mnt[1]){
                break
            } 
            else{
                vec_long[u] <- long_mnt[u] 
            }
         }
        aprendre <- seq(from=1,to=length(lat_mnt),by=length(vec_long)) 
        vec_lat <- lat_mnt[aprendre]
         
         he_mnt <- matrix(he_mnt,nrow=length(vec_long),byrow=FALSE)
         mat <- he_mnt
         if(vec_long[order(vec_long)[1]]!=vec_long[1]){
             vec_long <-  vec_long[order(vec_long)]
             longueur <- length(vec_long)
             for(u in 1:longueur){
                mat[u,] <- he_mnt[longueur-u+1,]
             }
         }
         he_mnt <- mat
         mat <- he_mnt
         if(vec_lat[order(vec_lat)[1]]!=vec_lat[1]){
             vec_lat <-  vec_lat[order(vec_lat)]
             longueur <- length(vec_lat)
             for(u in 1:longueur){
                  mat[,u] <- he_mnt[,longueur-u+1]
             }
         }
         he_mnt <- mat
        
        
        alti <- bicubic(x=vec_long,y=vec_lat,z=he_mnt,x0=long,y0=lat)$z

   }
  return(alti)
}




zoom <- function(mnt,X,Y,point,taille){ # Function which permits to zoom excluding members of a matrix farther than a fixed tolerance
  matrice_distance <- sqrt((point[1]-X)^2+(point[2]-Y)^2)
  mnt_local <- mnt[matrice_distance<taille]
  X_local <- X[matrice_distance<taille]
  Y_local <- Y[matrice_distance<taille]
  return(list(mnt_local,X_local,Y_local))
}





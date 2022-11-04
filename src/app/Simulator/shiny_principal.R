library(splus2R) # Loading of the splus2R library
library(pixmap)

source('pos2geo.R')
source('geo2pos.R')
source('projection.R')
source('dist_plan.R')
source('equ_plan.R')
source('reflexion_algo_gleason.R')
source('reflexion_algo_helm.R')
source('equ_ellipsoide.R')
source("coord_pts_tangente.R")
source("sp3.R")
source("read_parameter.R")
# STATIC CASE ONLY
principal <- function(input) {
  parameter.update.from.input(input)
  sp3.file = param["sp3_file_path"]
  sp3 = SP3.read(sp3.file)[1:1000, ]
  pos_E_xyz = cbind(as.numeric(sp3[, 2]), as.numeric(sp3[, 3]), as.numeric(sp3[, 4]))
  pos_E_xyz = asplit(pos_E_xyz, MARGIN = 1)
  
  # Calculation of the ellipsoid similar to the WGS84 ellipsoid, going through the ground projection of the receiver (= receiver with a zero altitude)
  pos_R_geo = as.numeric(c(lon = param["receiver_longitude"], lat = param["receiver_latitude"], he = param["receiver_height"]))
  pos_R_xyz = geo2pos(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_R_geo)
  ellipsoide_sol <-
    equ_ellipsoide(as.numeric(param["semi_major_WGS84"]),
                   as.numeric(param["semi_minor_WGS84"]),
                   pos_R_xyz,
                   as.numeric(param["receiver_eheight"])) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the ground projection of the receiver
  R_0_sol <- ellipsoide_sol[1]
  R_p_sol <- ellipsoide_sol[2]
  
  # Calculation of the ellipsoid similar to the WGS84, going through the receiver
  
  ellipsoide_recepteur <-
    equ_ellipsoide(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_R_xyz, 0) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the receiver
  R_0_recepteur <- ellipsoide_recepteur[1]
  R_p_recepteur <- ellipsoide_recepteur[2]
  
  # Calculation of the WGS84 cartesian coordinates of the ground projection of the receiver (= receiver with a zero altitude)
  pos_Rsol_xyz <- proj_surf(R_0_sol, R_p_sol, pos_R_xyz)
  pos_Rsol_geo <-
    proj_surf(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_Rsol_xyz)
  pos_Rellips_xyz <-
    proj_surf(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_R_xyz)
  pos_Rellips_geo <-
    pos2geo(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_Rellips_xyz)
  
  equ_plan_local <-
    equ_plan(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_R_xyz)
  
  x <- pos_Rsol_xyz[1]
  y <- pos_Rsol_xyz[2]
  x_centre <-
    x[1] # Center of the graph= position of the receiver point.
  y_centre <- y[1]
  coord_pts_tangente <-
    coord_pts_tangente(
      R_0_sol,
      R_p_sol,
      pos_R_xyz,
      pos_Rsol_xyz,
      pos2geo(R_0_sol, R_p_sol, pos_Rsol_xyz),
      pos2geo(R_0_sol, R_p_sol, pos_R_xyz)
    ) # coordinates of the two points of tangence of the radius going through the reception point, with the meridian ellipse
  dist_plani <- sqrt(pos_R_xyz[1] ^ 2 + pos_R_xyz[2] ^ 2)
  lim_max_1 <-
    abs(dist_plani - coord_pts_tangente[1]) # distance between the first tangence point with the receiver
  lim_max_2 <- abs(dist_plani - coord_pts_tangente[3])
  lim_max <-
    max(lim_max_1, lim_max_2) # the maximal distance between the reflection points with the receiver is the maximum of the two distance between the tangence points with the receiver
  lim_max_projete_1 <-
    sqrt((coord_pts_tangente[1] - dist_plani) ^ 2 + (coord_pts_tangente[2] -
                                                       pos_R_xyz[3]) ^ 2) # distance between the first tangence point with the receiver, projected in the local plane
  lim_max_projete_2 <-
    sqrt((coord_pts_tangente[3] - dist_plani) ^ 2 + (coord_pts_tangente[4] -
                                                       pos_R_xyz[3]) ^ 2)
  lim_max_projete <-
    max(lim_max_projete_1, lim_max_projete_2)# Maximal distance projected in the local plane
  if (is.na(lim_max)) {
    lim_max <- 10000
  }
  x_max <-  x_centre + lim_max / 100
  x_min <- x_centre - lim_max / 100
  y_max <-  y_centre + lim_max / 100
  y_min <- y_centre - lim_max / 100
  
  if (!input$checkbox_DEM) {
    # if we do not take the DEM into account
    
    # Coordinates in the "local" ellipsoid
    pos_Rsol_xyz <- proj_surf(R_0_sol, R_p_sol, pos_R_xyz)
    pos_Rsol_geo <- pos2geo(R_0_sol, R_p_sol, pos_Rsol_xyz)
    
    # Parameters of the osculating sphere
    R_gaussian <-
      R_0_sol ^ 2 * R_p_sol / ((R_0_sol * cos(pos_Rsol_geo[2] * pi / 180)) ^
                                 2 + (R_p_sol * sin(pos_Rsol_geo[2] * pi / 180)) ^ 2)
    pos_centre_geo <-
      c(pos_Rsol_geo[1], pos_Rsol_geo[2],-R_gaussian)
    pos_centre_xyz <- geo2pos(R_0_sol, R_p_sol, pos_centre_geo)
    
    # We calculate the satellite elevation
    
    vec_normal <- pos_Rsol_xyz - pos_centre_xyz
    # equation of the normal plane going through the receiver
    a <- vec_normal[1]
    b <- vec_normal[2]
    c1 <- vec_normal[3]
    d <- -(a * pos_R_xyz[1] + b * pos_R_xyz[2] + c1 * pos_R_xyz[3])
    distance_plan <-
      abs(unlist(lapply(
        pos_E_xyz, dist_plan, equ_plan = c(a, b, c1, d)
      )))
    sin_elevation_satellite <-
      abs(distance_plan / unlist(lapply(pos_E_xyz, function(x) {
        vecnorm(x - c(pos_R_xyz))
      }))) #vecnorm(pos_E_xyz - pos_R_xyz))
    sin_elevation_satellite[sin_elevation_satellite >= 1] = -1
    elevation_satellite <- asin(sin_elevation_satellite)
    
    # equation of the local plane
    d <-
      -(a * pos_Rsol_xyz[1] + b * pos_Rsol_xyz[2] + c1 * pos_Rsol_xyz[3])
    distance_plan <-
      abs(unlist(lapply(pos_E_xyz, FUN = dist_plan, equ_plan = c(a, b, c1, d))))
    # normal distance from the satellite to the local plane
    if (param["algorithm"] == 'ellipsoid') {
      # If the user wants to use the ellipsoid algorithm
      source('calcul_normale.R')
      # si on veut tester l'algo en mettant les caractoristiques de la spheree
      #               pos_Ei_xyz <- pos_E_xyz -pos_centre_xyz
      #                  pos_Ri_xyz <- pos_R_xyz -pos_centre_xyz
      #                  pos_S <- reflexion_algo_gleason(R_gaussian,R_gaussian,pos_Ei_xyz,pos_Ri_xyz,as.numeric(param[19,2]))
      #                  pos_S_xyz <- c(pos_S[[2]][1],pos_S[[3]][1],pos_S[[4]][1])+pos_centre_xyz
      pos_S <-
        reflexion_algo_gleason(R_0_sol,
                               R_p_sol,
                               pos_E_xyz,
                               pos_R_xyz,
                               as.numeric(param["tolerance"]))
      pos_S_xyz <- c(pos_S[[2]][1], pos_S[[3]][1], pos_S[[4]][1])
      pos_S_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
    }
    
    else if (param["algorithm"] == 'plane') {
      source('reflexion_algo_plan.R')
      pos_R_geo <- pos2geo(R_0_sol, R_p_sol, pos_R_xyz)
      pos_Rfictif_geo <-
        c(pos_R_geo[1],
          pos_R_geo[2],
          pos_R_geo[3] - 2 * as.numeric(param["receiver_height"]))
      pos_Rfictif_xyz <- geo2pos(R_0_sol, R_p_sol, pos_Rfictif_geo)
      pos_Esol_xyz <-
        lapply(pos_E_xyz, proj_plan, equ_plan = c(a, b, c1, d)) #, pos_E_xyz)
      pos_S <-
        mcmapply(
          reflexion_algo_plan,
          pos_E_xyz = pos_E_xyz,
          pos_Esol_xyz = pos_Esol_xyz,
          distance_plan = distance_plan,
          pos_R_xyz,
          MoreArgs = list(
            pos_Rsol_xyz,
            pos_Rfictif_xyz,
            as.numeric(param["receiver_height"]),
            a,
            b,
            c1,
            d
          ), mc.cores = 4
        ) # We calculate the reflection point position
      pos_S_xyz <- t(pos_S)[, c(2, 3, 4)]
      pos_S_xyz = as.data.frame(pos_S_xyz)
      pos_S_geo <- apply(pos_S_xyz, 1, pos2geo, R_0 = as.numeric(param["semi_major_WGS84"]),
                           R_p = as.numeric(param["semi_minor_WGS84"]))
      closeAllConnections()
      pos_S_geo = t(pos_S_geo)
    }
    else if (param["algorithm"] == 'sphere') {
      # If the user wants to use the sphere algorithm
      pos_S <- mcmapply(reflexion_algo_helm,
          R_gaussian,
          pos_centre_xyz,
          pos_E_xyz,
          pos_R_xyz,
          elevation_satellite,
          as.numeric(param["receiver_height"]),
          param["tolerance"]
        )
      
      pos_S_geo <-
        apply(pos_S_xyz,
              1,
              pos2geo,
              R_0 = as.numeric(param["semi_major_WGS84"]),
              R_p = as.numeric(param["semi_minor_WGS84"]))
      R_sphere <- as.numeric(pos_S[7])
    }
  }
  return(tibble::add_column(data.frame(pos_S_geo), sp3[,1], ID = sp3[,"ID"]))
}

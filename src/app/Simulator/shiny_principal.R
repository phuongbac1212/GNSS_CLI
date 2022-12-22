library(splus2R) # Loading of the splus2R library
library(pixmap)
library(DescTools)

source('pos2geo.R')
source('geo2pos.R')
source('projection.R')
source('dist_plan.R')
source('equ_plan.R')
source('reflexion_algo_gleason.R')
source('reflexion_algo_helm.R')
source('equ_ellipsoide.R')
source("coord_pts_tangente.R")
source("calculate_azimuth_elevation.R")
source("sp3.R")
source("surf_Fresnel.r")
source("read_parameter.R")
source("visibilite_sat.r")
source("determination_coord.R")
# require(xlsx)
# STATIC CASE ONLY
principal <- function(input) {
  print("start cal")
  # parameter.update.from.input(input)
  sp3.file = param["sp3_file_path"]
  sp3 = data.frame(SP3.read(sp3.file))
  #sp3$time = lubridate::as_datetime(as.numeric(sp3$TimeStamp), origin = "1980-1-6")
  sp3 = filter(sp3, ID == "PG01")[1:80,]
  # sp3 = filter(sp3, ID == "PG01")
  pos_R_geo = as.numeric(c(lon = param["receiver_longitude"], lat = param["receiver_latitude"], he = param["receiver_height"]))
  pos_R_xyz = geo2pos(as.numeric(param["semi_major_WGS84"]), as.numeric(param["semi_minor_WGS84"]), pos_R_geo)
  
  # Calculation of the ellipsoid similar to the WGS84 ellipsoid, going through the ground projection of the receiver (= receiver with a zero altitude)
  
  
  ellipsoide_sol <-
    equ_ellipsoide(
      R_0_WGS84 = as.numeric(param["semi_major_WGS84"]),
      R_p_WGS84 = as.numeric(param["semi_minor_WGS84"]),
      pos_xyz = pos_R_xyz,
      hauteur = as.numeric(param["receiver_height"])
    ) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the ground projection of the receiver
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
  
  
  pos_E_xyz = cbind(as.numeric(sp3[, 2]), as.numeric(sp3[, 3]), as.numeric(sp3[, 4]))
  pos_E_xyz = asplit(pos_E_xyz, MARGIN = 1)
  
  vis = mapply(sat_visible, pos_E_xyz, MoreArgs = list(
    pos_R_xyz = pos_R_xyz,
    R_0_WGS84 = as.numeric(param["semi_major_WGS84"]), 
    R_p_WGS84 = as.numeric(param["semi_minor_WGS84"]),
    masque_vert_min = as.numeric(param["low_vertical_mask"]),
    masque_vert_max = as.numeric(param["up_vertical_mask"]),
    masque_hor_min = as.numeric(param["low_horizontal_mask"]),
    masque_hor_max = as.numeric(param["up_horizontal_mask"]),
    equ_plan_local = equ_plan_local
  )) %>% t() %>% data.frame()
  names(vis) = c("visible", "azimuth", "elevation")
  
  sp3 = tibble::add_column(sp3, vis) %>% filter(visible == T)
  
  pos_E_xyz = cbind(as.numeric(sp3[, 2]), as.numeric(sp3[, 3]), as.numeric(sp3[, 4]))
  pos_E_xyz = asplit(pos_E_xyz, MARGIN = 1)
  elevation_satellite=NA
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
      abs(unlist(lapply(
        pos_E_xyz, FUN = dist_plan, equ_plan = c(a, b, c1, d)
      )))
    # normal distance from the satellite to the local plane
    if (param["algorithm"] == 'ellipsoid') {
      print("ellipsoid")
      # If the user wants to use the ellipsoid algorithm
      source('calcul_normale.R')
      # si on veut tester l'algo en mettant les caractoristiques de la spheree
      #               pos_Ei_xyz <- pos_E_xyz -pos_centre_xyz
      #                  pos_Ri_xyz <- pos_R_xyz -pos_centre_xyz
      #                  pos_S <- reflexion_algo_gleason(R_gaussian,R_gaussian,pos_Ei_xyz,pos_Ri_xyz,as.numeric(param[19,2]))
      #                  pos_S_xyz <- c(pos_S[[2]][1],pos_S[[3]][1],pos_S[[4]][1])+pos_centre_xyz
      pos_S <-
        mclapply(X = pos_E_xyz, 
               FUN= reflexion_algo_gleason,
               R_0 = R_0_sol,
               R_p = R_p_sol,
               pos_R_xyz= pos_R_xyz,
               tolerance_gleason = as.numeric(param["tolerance"]), mc.cores = 4) %>% rbindlist()
      pos_S = filter(pos_S, V1 == T)
      pos_S_xyz <- pos_S[, 2:4]
      pos_S_geo <-
        apply(pos_S_xyz,
              1,
              pos2geo,
              R_0 = as.numeric(param["semi_major_WGS84"]),
              R_p = as.numeric(param["semi_minor_WGS84"])) %>% t()
      #write.xlsx(pos_S_geo, "a.xlsx", sheetName = "ellipsoid", append = T)
    }
    
    else if (param["algorithm"] == 'plane') {
      source('reflexion_algo_plan.R')
      print("Plane")
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
          ),
          mc.cores = 4
        ) # We calculate the reflection point position
      closeAllConnections()
      
      pos_S_xyz <- t(pos_S)[, c(2, 3, 4)]
      pos_S_xyz = as.data.frame(pos_S_xyz)
      pos_S_geo <-
        apply(pos_S_xyz,
              1,
              pos2geo,
              R_0 = as.numeric(param["semi_major_WGS84"]),
              R_p = as.numeric(param["semi_minor_WGS84"]))
      pos_S_geo = t(pos_S_geo)
      #write.xlsx(pos_S_geo, "a.xlsx", sheetName = "plane", append = T)
    }
    else if (param["algorithm"] == 'sphere') {
      print("sphere")
      # If the user wants to use the sphere algorithm
      pos_S_xyz <- mcmapply(
        reflexion_algo_helm,
        pos_E_xyz = pos_E_xyz,
        elevation = elevation_satellite,
        MoreArgs = list(
          R_T = R_gaussian,
          pos_centre_xyz = pos_centre_xyz,
          pos_R_xyz = pos_R_xyz,
          hauteur = as.numeric(param["receiver_height"]),
          tolerance_helm = param["tolerance"]
        )
      )
      pos_S_xyz = data.frame(t(pos_S_xyz))[, c(2, 3, 4)]
      pos_S_geo <-
        apply(pos_S_xyz,
              1,
              pos2geo,
              R_0 = as.numeric(param["semi_major_WGS84"]),
              R_p = as.numeric(param["semi_minor_WGS84"]))
      pos_S_geo = t(pos_S_geo)
      #R_sphere <- as.numeric(pos_S[7])
     # write.xlsx(pos_S_geo, "a.xlsx", sheetName = "sphere", append = T)
    }
  }
  wavelength = 299792458 / switch(param["wavelength"],
                                  "L1" = 1575420000,
                                  "L2" = 1227600000,
                                  "L5" = 1176000000)
  result_fresnel <- mapply(
    surf_Fresnel,
    lambda = wavelength,
    elev = elevation_satellite,
    MoreArgs = list(hauteur = as.numeric(param["receiver_height"]))
  )
  
  res = tibble::add_column(data.frame(t(result_fresnel)) ,
                           data.frame(pos_S_xyz),
                           sp3)
  names(res)[1:6] = c("a", "b", "c", "XS", "YS", "ZS")
  res$color = unlist(apply(res,1,FUN = function(x) {
    x = substr(x[["ID"]], 2,2)
    return(switch(x,
                  G = "red",
                  R = "blue",
                  E = "yellow",
                  C = "green"))
  }))
  
  gcp = apply(res, 1, function(x) {
    pos_E_xyz = as.numeric(c(x[["X"]],x[["Y"]], x[["Z"]]))
    pos_S_xyz = as.numeric(c(x[["XS"]],x[["YS"]], x[["ZS"]]))
    equ_plan_local_S <- equ_plan(pos_xyz =  pos_S_xyz)
    plan_pos_R_xyz <- proj_plan(equ_plan_local_S, pos_R_xyz)
    plan_pos_E_xyz <- proj_plan(equ_plan_local_S, pos_E_xyz)
    plan_pos_S_xyz <- proj_plan(equ_plan_local_S, pos_S_xyz)
    vec_droite <-
      c(
        plan_pos_R_xyz[1] - plan_pos_S_xyz[1],
        plan_pos_R_xyz[2] - plan_pos_S_xyz[2],
        plan_pos_R_xyz[3] - plan_pos_S_xyz[3]
      )
    k <- result_fresnel[1] / vecnorm(vec_droite)
    coord_A <-
      c(
        vec_droite[1] * k + plan_pos_S_xyz[1],
        vec_droite[2] * k + plan_pos_S_xyz[2],
        vec_droite[3] * k + plan_pos_S_xyz[3]
      )
    k <- -result_fresnel[1] / vecnorm(vec_droite)
    coord_AA <-
      c(
        vec_droite[1] * k + plan_pos_S_xyz[1],
        vec_droite[2] * k + plan_pos_S_xyz[2],
        vec_droite[3] * k + plan_pos_S_xyz[3]
      )
    
    # We determine the coordinates of B
    
    Xa <- coord_A[1]
    Ya <- coord_A[2]
    Za <- coord_A[3]
    u  <- equ_plan_local_S[1]
    v <- equ_plan_local_S[2]
    w <- equ_plan_local_S[3]
    r <- equ_plan_local_S[4]
    Xo <- plan_pos_S_xyz[1]
    Yo <- plan_pos_S_xyz[2]
    Zo <- plan_pos_S_xyz[3]
    b <- result_fresnel[2]
    
    result <-
      determination_coord(Xa, Ya, Za, u, v, w, r, Xo, Yo, Zo, b)
    coord_B1 <- c(result[1], result[2], result[3])
    coord_B2 <- c(result[4], result[5], result[6])
    
    vec_S_B1 <-
      c(
        coord_B1[1] - plan_pos_S_xyz[1],
        coord_B1[2] - plan_pos_S_xyz[2],
        coord_B1[3] - plan_pos_S_xyz[3]
      )
    vec_S_B2 <-
      c(
        coord_B2[1] - plan_pos_S_xyz[1],
        coord_B2[2] - plan_pos_S_xyz[2],
        coord_B2[3] - plan_pos_S_xyz[3]
      )
    vec_S_R <-
      c(
        plan_pos_R_xyz[1] - plan_pos_R_xyz[1],
        plan_pos_R_xyz[2] - pos_S_xyz[2],
        plan_pos_R_xyz[3] - plan_pos_S_xyz[3]
      )
    
    det_1 <-  vec_S_R[1] * vec_S_B1[2] - vec_S_R[2] * vec_S_B1[1]
    det_2 <-  vec_S_R[1] * vec_S_B2[2] - vec_S_R[2] * vec_S_B2[1]
    
    sin_1 <- det_1 / vecnorm(vec_S_R) / vecnorm(vec_S_B1)
    sin_2 <- det_2 / vecnorm(vec_S_R) / vecnorm(vec_S_B2)
    
    if (sin_1 > 0) {
      coord_B <- coord_B1
    }
    else{
      coord_B <- coord_B2
    }
    
    ellipsoide_b <-
      equ_ellipsoide(pos_xyz = coord_AA, hauteur = result_fresnel[2]) # We get R_0 and R_p (semi-major and -minor axis of the ellipsoid parallel to the WGS84 and going through the ground point
    R_0_b <- ellipsoide_b[1]
    R_p_b <- ellipsoide_b[2]
    pos_AA_BB_xyz <- proj_surf(R_0_b, R_p_b, coord_AA)
    coord_AA_BB <- pos_AA_BB_xyz
    names(coord_A)  = c("XA" , "YA", "ZA");
    names(coord_B)  = c("XB" , "YB", "ZB");
    names(coord_AA_BB)  = c("XAABB" , "YAABB", "ZAABB");
    return(c(coord_A,
           coord_B,
           coord_AA_BB))
  })
  gcp = data.frame(t(gcp))
  res = tibble::add_column(res, gcp)
  print("end calculation")
  return(res)
}

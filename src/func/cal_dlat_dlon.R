source("src/func/equ_ellipsoide.R")
source("src/func/calcul_normale.R")
source("src/func/dist_plan.R")
source("src/func/projection.R")
source("src/func/determination_coord.R")
source("src/func/surf_Fresnel.R")
source("src/func/equ_plan.R")
source("src/func/pos2geo.R")
source("src/func/geo2pos.R")
cal.dlat.dlon <-
  function(pos_E_xyz,
           pos_R_xyz,
           elevation,
           azimuth,
           wavelength) {
    longueur_onde = wavelength
    # pos = pos2geo(pos_xyz = c(receiver.coor['x'], receiver.coor["y"], receiver.coor["z"]))
    # dist_rp = h0 / sin(as.numeric(elevation))
    # dx = dist_rp * sin(as.numeric(azimuth))
    # dy = dist_rp * cos(as.numeric(azimuth))
    # dLon = dx / (111320 * cos(pos[2]))
    # dLat = dy / 110540
    # lat_rp = pos[2] + dLat
    # lon_rp = pos[1] - dLon
    # h_rp = pos[3] - h0
    ellipsoide_sol <-
      equ_ellipsoide(pos_xyz = pos_R_xyz, hauteur = h0) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the ground projection of the receiver
    R_0_sol <- ellipsoide_sol[1]
    R_p_sol <- ellipsoide_sol[2]
    # if we do not take the DEM into account
    # Coordinates in the "local" ellipsoid
    pos_Rsol_xyz <- proj_surf(R_0 = R_0_sol, R_p = R_p_sol, pos_R_xyz)
    pos_Rsol_geo <- pos2geo(R_0 = R_0_sol, R_p = R_p_sol, pos_Rsol_xyz)
    
    # Parameters of the osculating sphere
    R_gaussian <-
      R_0_sol ^ 2 * R_p_sol / ((R_0_sol * cos(pos_Rsol_geo[2] * pi / 180)) ^
                                 2 + (R_p_sol * sin(pos_Rsol_geo[2] * pi / 180)) ^ 2)
    pos_centre_geo <-
      cbind(pos_Rsol_geo[1], pos_Rsol_geo[2], -R_gaussian)
    pos_centre_xyz <-
      geo2pos(R_0 = R_0_sol, R_p = R_p_sol, pos_centre_geo)
    
    # We calculate the satellite elevation
    
    vec_normal <- pos_Rsol_xyz - pos_centre_xyz
    # equation of the normal plane going through the receiver
    a <- vec_normal[1]
    b <- vec_normal[2]
    c1 <- vec_normal[3]
    d <- -(a * pos_R_xyz[1] + b * pos_R_xyz[2] + c1 * pos_R_xyz[3])
    distance_plan <-
      abs(dist_plan(c(a, b, c1, d), pos_E_xyz)) # We calculate the distance between the satellite and its projection w.r.t the plane going through the receiver, tangent to the ellipsoid going throught the receiver
    sin_elevation_satellite <-
      abs(distance_plan / vecnorm(pos_E_xyz - pos_R_xyz))
    if (sin_elevation_satellite >= 1) {
      sin_elevation_satellite <- 1
    }
    elevation_satellite <- asin(sin_elevation_satellite)
    
    # equation of the local plane
    d <-
      -(a * pos_Rsol_xyz[1] + b * pos_Rsol_xyz[2] + c1 * pos_Rsol_xyz[3])
    distance_plan <-
      abs(dist_plan(c(a, b, c1, d), pos_E_xyz)) # normal distance from the satellite to the local plane
    source('src/func/reflexion_algo_plan.R')
    pos_R_geo <- pos2geo(pos_R_xyz, R_0_sol, R_p_sol)
    pos_Rfictif_geo <-
      cbind(pos_R_geo[1], pos_R_geo[2], pos_R_geo[3] - 2 * h0)
    pos_Rfictif_xyz <- geo2pos(pos_Rfictif_geo, R_0_sol, R_p_sol)
    pos_Esol_xyz <- proj_plan(c(a, b, c1, d), pos_E_xyz)
    pos_S <-
      reflexion_algo_plan(
        pos_E_xyz,
        pos_Esol_xyz,
        pos_R_xyz,
        pos_Rsol_xyz,
        pos_Rfictif_xyz,
        h0,
        distance_plan,
        a,
        b,
        c1,
        d
      ) # We calculate the reflection point position
    pos_S_xyz <- c(pos_S[[2]][1], pos_S[[3]][1], pos_S[[4]][1])
    pos_S_geo <- pos2geo(pos_S_xyz)
    
    
    # cal fresnel zone
    # Fresnel surface
    result_fresnel <-
      surf_Fresnel(longueur_onde, h0, elevation_satellite)
    
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
    
    result <- determination_coord(Xa, Ya, Za, u, v, w, r, Xo, Yo, Zo, b)
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
    res = c(pos_S_geo, result_fresnel, coord_A, coord_B, coord_AA_BB)
    names(res) = c(
      "Lon_simu",
      "Lat_simu",
      "Height_simu",
      "Fz_a",
      "Fz_b",
      "Fz_c",
      "A_x",
      "A_y",
      "A_z",
      "B_x",
      "B_y",
      "B_z",
      "AB_x",
      "AB_y",
      "AB_z"
    )
    return(res)
  }

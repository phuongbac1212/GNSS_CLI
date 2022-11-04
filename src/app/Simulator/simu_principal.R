##################################################################
#
# Main program of the simulator
#
# launch: source("simu_principal.R")
#
# 2nd semester 2012
#
# Author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

library(splus2R) # Loading of the splus2R library
library(pixmap)

# Loading of the differents scripts or functions which will be used
#

source('pos2geo.R')
source('geo2pos.r')
source('projection.r')
source('dist_plan.r')
source('equ_plan.r')
source('reflexion_algo_gleason.R')
source('reflexion_algo_helm.r')
source('equ_ellipsoide.r')
source('coord_pts_tangente.r')
source('tangente_ellipse.r')
source('txt2kml.r')
source('visibilite_sat.r')
source('mnt_fonctions.r')
source('reflexion_mnt.r')
source('conversion_alti_he.r')
source('surf_Fresnel.r')
source('determination_coord.r')
require(tcltk)

options(digits = 22)   # To get enough digits for calculations

principal <-
  function(nom_fichier_config,
           mon_data,
           nombre_satellites,
           nb_ligne) {
    T1 <- Sys.time()
    system(
      'wine simu_del_all.bat',
      show.output.on.console = FALSE,
      ignore.stderr = TRUE,
      ignore.stdout = TRUE
    ) # We delete all the old files of the output folder
    nom_fichier_pts_reflexion_glonass <-  'specular_points_glonass.txt'
    nom_fichier_sat_glonass <-  "satellites_glonass.txt"
    nom_fichier_pts_reflexion_gps <-  'specular_points_gps.txt'
    nom_fichier_sat_gps <- "satellites_gps.txt"
    nom_fichier_pts_reflexion_galileo <-  'specular_points_galileo.txt'
    nom_fichier_sat_galileo <-    "satellites_galileo.txt"
    
    write.table(
      'Longitude (o) - Latitude (o) - Ellipdoidal height (o) - Date - PRN - Azimuth angle (o) - Elevation angle (o) - Incr',
      gsub(' ', '', paste(path_output, nom_fichier_sat_glonass)),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    write.table(
      'Longitude (o) - Latitude (o) - Ellipdoidal height (o) - Date - PRN - Azimuth angle (o) - Elevation angle (o) - Incr',
      gsub(' ', '', paste(path_output, nom_fichier_sat_gps)),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    write.table(
      'Longitude (o) - Latitude (o) - Ellipdoidal height (o) - Date - PRN - Azimuth angle (o) - Elevation angle (o) - Incr',
      gsub(' ', '', paste(path_output, nom_fichier_sat_galileo)),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    write.table(
      'Longitude (o) - Latitude (o) - Ellipdoidal height (o) - Date - Incr - PRN - Azimuth angle (o) - Elevation angle (o) - Semi major axis First Fresnel area (m) - Semi minor axis First Fresnel area (m) - Area of the first Fresnel Area (m2)',
      gsub(
        ' ',
        '',
        paste(path_output, nom_fichier_pts_reflexion_glonass)
      ),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    write.table(
      'Longitude (o) - Latitude (o) - Ellipdoidal height (o) - Date - Incr - PRN - Azimuth angle (o) - Elevation angle (o) - Semi major axis First Fresnel area (m) - Semi minor axis First Fresnel area (m) - Area of the first Fresnel Area (m2)',
      gsub(' ', '', paste(
        path_output, nom_fichier_pts_reflexion_gps
      )),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    write.table(
      'Longitude (o) - Latitude (o) - Ellipdoidal height (o) - Date - Incr - PRN - Azimuth angle (o) - Elevation angle (o) - Semi major axis First Fresnel area (m) - Semi minor axis First Fresnel area (m) - Area of the first Fresnel Area (m2)',
      gsub(
        ' ',
        '',
        paste(path_output, nom_fichier_pts_reflexion_galileo)
      ),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    write.table(
      'Longitude (o) - Latitude (o) - Ellipsoidal height (m) - vertical mask (o) - horizontal mask (o) - Date',
      gsub(' ', '', paste(path_output, "receiver.txt")),
      append = TRUE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    
    
    param <-
      as.matrix(read.table(
        nom_fichier_config,
        header = FALSE,
        sep = "",
        na.strings = "NA",
        fill = TRUE
      ))
    
    nom_fichier <-
      as.character(param[1, 2]) # nom_fichier : name of the ephemeris file. Pattern:  'path/nom.xxx'
    static_receiver <<- TRUE
    if (grep('.pos', c(nom_fichier, '.pos'))[1] == 1) {
      static_receiver <<- FALSE
    }
    nom_mnt <-
      as.character(param[2, 2]) # nom_mnt : name of the DEM file. Pattern: 'path/nom.xxx'
    heure_deb <-
      as.numeric(param[3, 2]) # heure_deb, min_deb : hour of the beginning of the simulation
    min_deb <- as.numeric(param[4, 2])
    heure_fin <-
      as.numeric(param[5, 2]) # heure_fin, min_fin : hour of the end of the simulation
    min_fin <- as.numeric(param[6, 2])
    hauteur_recepteur_init <-
      as.numeric(param[7, 2]) # hauteur_recepteur : the receiver height with respect to the ground surface
    hauteur_recepteur <- hauteur_recepteur_init
    R_0_WGS84 <-
      as.numeric(param[21, 2]) # R_0_WGS84 : semi-major axis of the WGS84 ellipsoid (m)
    f_WGS84 <-
      1 / as.numeric(param[22, 2]) # f_WGS84 : flatness of the WGS84 ellipsoid
    R_p_WGS84 <-
      R_0_WGS84 * (1 - f_WGS84) # Semi-minor axis of the WGS84 ellipsoid (m)
    pos_R_long <- as.numeric(param[8, 2])
    pos_R_lat  <- as.numeric(param[9, 2])
    pos_R_he <- as.numeric(param[10, 2])
    pos_R_xyz <-
      geo2pos(R_0_WGS84, R_p_WGS84, c(pos_R_long, pos_R_lat, pos_R_he)) # pos_R_xyz : cartesian coordinates of the receiver in the WGS84 system
    masque_vert_min <- as.numeric(param[12, 2]) # vertical mask
    masque_vert_max <- as.numeric(param[13, 2])
    masque_hor_min <- as.numeric(param[14, 2]) # horizontal mask
    masque_hor_max <- as.numeric(param[15, 2])
    glonass  <-
      as.numeric(param[16, 2]) # do we intergrate the GLONASS constellation o
    algorithme <-
      as.numeric(param[17, 2]) # Which algorithm do we take to calculate the reflection o
    facteur_K <- as.numeric(param[18, 2]) # Value of the K factor
    tolerance_calcul_ref_gleason <-
      as.numeric(param[19, 2]) # tolerance used in the reflection calculation (ellipsoid algorithm)
    tolerance_calcul_ref_helm <-
      as.numeric(param[20, 2]) # tolerance used in the reflection calculation (sphere aglorithm)
    date_simu <- param[23, 2] # Date of the simulation
    longueur_onde <- param[34, 2] # wave length
    dessiner_fresnel <-
      param[35, 2] # determination of the first Fresnel surface o
    afficher_output <-
      param[36, 2] # display of the progression (take more time) o
    corriger_tropo <-
      param[37, 2] # correction of the tropospheric error o
    nom_fichier_amfx_sat <-
      param[38, 2] # name of the AMF file calculated from the ground until the end of the tropo (6450 km)
    nom_fichier_amfx_sta <-
      param[39, 2]#  name of the AMF file calculated from the ground until the receiver
    nom_station_amfx <-
      param[40, 2] # name of the station in the AMF files
    if (dessiner_fresnel == "Yes") {
      dessiner_fresnel <- TRUE
    }
    else{
      dessiner_fresnel <- FALSE
    }
    if (afficher_output == 'Yes') {
      affichage_output <- TRUE
    }
    else{
      affichage_output <- FALSE
    }
    if (corriger_tropo == 'yes') {
      correction_tropo <- TRUE
    }
    else{
      correction_tropo <- FALSE
    }
    # mon_data : tab containing the parameters read in the epheremeris file ==> cf. lect_fichier.R
    # nombre_satellites : number of satellites whose coordinates are known (in the ephemeris file)
    # nb_ligne : number of lines read in the ephemeris file
    
    #
    # Parameters describing the WGS84 ellipsoid
    #
    
    
    e_WGS84 <-
      sqrt((R_0_WGS84 ^ 2 - R_p_WGS84 ^ 2) / R_0_WGS84 ^ 2)  # excentricity of the WGS84 ellipsoid
    
    #
    # Determination of the Work Directory and the different paths
    #
    
    #wd <- getwd()
    #setwd('..')
    #path_output <- gsub(' ','',paste(getwd(),'/Output/Simulation/')) # path of the Output folder
    #path_input_pos_sat <- gsub(' ','',paste(getwd(),'/Input/Satellites_position/')) # path of the Input folder
    #path_input_amfx <- gsub(' ','',paste(getwd(),'/Input/AMFX/')) # Path of the AMFX folder
    #path_input_mnt <- gsub(' ','',paste(getwd(),'/Input/DEM/')) # Path of the DEM folder
    #setwd(gsub(' ','',paste(getwd(),'/Programs'))) # We put back the Work Directory in the Programs folder
    
    #
    # Preparation of the output
    #
    
    write.table(
      "",
      file = gsub(' ', '', paste(path_output, "output.txt")),
      append = FALSE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )      # Efface les donnoes de output.txt
    ecrit <-
      c(
        ' ********************************************************************************\n',
        ' GRSS - Time-varying positions of the specular reflection points \n',
        ' ********************************************************************************\n\n'
      )
    write.table(
      ecrit,
      file = gsub(' ', '', paste(path_output, "output.txt")),
      append = TRUE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    if (affichage_output) {
      write.table(
        ecrit,
        file = stdout(),
        append = TRUE,
        quote = FALSE,
        sep = " ",
        na = "-",
        row.names = FALSE,
        col.names = FALSE
      )
    }
    
    if (static_receiver) {
      #
      # Calculation of the ellipsoid similar to the WGS84 ellipsoid, going through the ground projection of the receiver (= receiver with a zero altitude)
      #
      
      ellipsoide_sol <-
        equ_ellipsoide(R_0_WGS84, R_p_WGS84, pos_R_xyz, hauteur_recepteur) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the ground projection of the receiver
      R_0_sol <- ellipsoide_sol[1]
      R_p_sol <- ellipsoide_sol[2]
      
      #
      # Calculation of the ellipsoid similar to the WGS84, going through the receiver
      #
      
      ellipsoide_recepteur <-
        equ_ellipsoide(R_0_WGS84, R_p_WGS84, pos_R_xyz, 0) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the receiver
      R_0_recepteur <- ellipsoide_recepteur[1]
      R_p_recepteur <- ellipsoide_recepteur[2]
      
      #
      # Calculation of the WGS84 cartesian coordinates of the ground projection of the receiver (= receiver with a zero altitude)
      #
      
      pos_Rsol_xyz <- proj_surf(R_0_sol, R_p_sol, pos_R_xyz)
      pos_Rsol_geo <-  proj_surf(R_0_WGS84, R_p_WGS84, pos_Rsol_xyz)
      pos_Rellips_xyz <- proj_surf(R_0_WGS84, R_p_WGS84, pos_R_xyz)
      pos_Rellips_geo <-
        pos2geo(R_0_WGS84, R_p_WGS84, pos_Rellips_xyz)
      
      
      equ_plan_local <- equ_plan(R_0_WGS84, R_p_WGS84, pos_R_xyz)
    }
    if (!static_receiver) {
      mon_data_pos <-
        scan(
          paste(path_input_pos_sat, nom_fichier, sep = ''),
          what = character(0),
          quiet = TRUE
        )
      coord_rec <-
        mon_data_pos[(which(mon_data_pos == 'ratio') + 1):length(mon_data_pos)]
      coord_rec <- matrix(coord_rec, ncol = 15, byrow = TRUE)
      ma_date <-
        as.POSIXct(paste(gsub('/', '-', coord_rec[, 1]), ' ', coord_rec[, 2], sep =
                           ''))
      lat_rec <- as.numeric(coord_rec[, 3])
      long_rec <- as.numeric(coord_rec[, 4])
      he_rec <- as.numeric(coord_rec[, 5])
      pos_R_geo_tot <- matrix(c(long_rec, lat_rec, he_rec), ncol = 3)
    }
    
    #
    # We inform the user that the ephemeris file has successfully been read
    #
    
    ecrit <-
      t(c(
        '',
        nb_ligne,
        'lines of the file',
        nom_fichier,
        'has been taken into account. \n'
      ))
    write.table(
      ecrit,
      file = gsub(' ', '', paste(path_output, "output.txt")),
      append = TRUE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    if (affichage_output) {
      write.table(
        ecrit,
        file = stdout(),
        append = TRUE,
        quote = FALSE,
        sep = " ",
        na = "-",
        row.names = FALSE,
        col.names = FALSE
      )
    }
    #
    # We display all the chosen parameters
    #
    
    ecrit <- gsub('_', ' ', param)
    write.table(
      ecrit,
      file = gsub(' ', '', paste(path_output, "output.txt")),
      append = TRUE,
      quote = FALSE,
      sep = " ",
      na = "-",
      row.names = FALSE,
      col.names = FALSE
    )
    if (affichage_output) {
      write.table(
        ecrit,
        file = stdout(),
        append = TRUE,
        quote = FALSE,
        sep = " ",
        na = "-",
        row.names = FALSE,
        col.names = FALSE
      )
    }
    
    
    #
    # initialisation of the variables
    #
    
    num_ligne <- 1
    compteur <- 0
    M <- 0
    L <- 0
    iter <- 1
    
    
    if (static_receiver) {
      x <- pos_Rsol_xyz[1]
      y <- pos_Rsol_xyz[2]
      
      #
      # We will plot the reflection point position in a graph whose axis are the WGS84 cartesian coordinates of the reflection points.
      #
      
      x_centre <-
        x[1] # Center of the graph= position of the receiver point.
      y_centre <- y[1]
      
      #
      # We want to determine xlim et ylim the min and max of the graph to display.
      # To do that, we know that the maximal reflection angle is 90o, which means that the wave is tangent to the ellipsoid
      # We will look for the coordinates of the points of tangence of the two radius going through the reception point with the meridian ellipse.
      # We will then calculate the distances between these two points and the reception point, which will correspond to the maximal distance between a reflection point and a given receiver
      # (This distance is indeed maximal for a radius tangent in the meridian plane. In the plane of the parallel, the section would be a circle and not any more an ellipse)
      #
      
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
    }
    
    #
    # Launchof the calculation
    #
    
    #
    # Preparation of the DEM if the user wants to take it into account
    #
    
    nom_mnt <- param[2, 2]
    mnt_rogne <- 0
    X_mnt <- 0
    Y_mnt <- 0
    if (nom_mnt != '-') {
      source('fenetre_info_mnt.r')
      prologue <-
        scan(
          paste(path_input_mnt, nom_mnt, sep = ""),
          what = character(0),
          sep = "",
          quiet = TRUE,
          nlines = 6
        )
      dimension <- as.numeric(prologue[2])
      x_BG_corner <- as.numeric(prologue[6])
      y_BG_corner <- as.numeric(prologue[8])
      cellsize <- as.numeric(prologue[10])
      x_BD_corner <- x_BG_corner + cellsize * dimension
      x_HD_corner <- x_BG_corner + cellsize * dimension
      y_HG_corner <- y_BG_corner + cellsize * dimension
      y_HD_corner <- y_BG_corner + cellsize * dimension
      cellsize <- as.numeric(prologue[10])
      mnt <-
        scan(
          paste(path_input_mnt, nom_mnt, sep = ""),
          what = numeric(0),
          sep = "",
          quiet = TRUE,
          skip = 6
        )
      mnt <- matrix(mnt, nrow = dimension , byrow = TRUE)
      
      
      tkdestroy(ttinfo)
      
      
      mnt_rogne <- 0
      ordre <- 0
      M_lim <- 0
      L_lim <- 0
      
      division <- param[30, 2]
      x_BG_check <-
        as.numeric(param[26, 2]) # abscissa of the DEM bottom left hand corner
      y_BG_check <-
        as.numeric(param[27, 2]) # ordinate of the DEM bottom left hand corner
      x_HD_check <-
        as.numeric(param[28, 2]) # abscissa of the DEM upper right hand corner
      y_HD_check <-
        as.numeric(param[29, 2]) # ordinate of the DEM upper right hand corner
      division <- as.numeric(param[30, 2])
      
      
      # We trim the DEM
      x_BG_lim <- trunc((x_BG_check - x_BG_corner) / cellsize)
      y_BG_lim <- trunc((y_BG_check - y_BG_corner) / cellsize)
      x_HD_lim <- trunc((x_HD_check - x_BG_corner) / cellsize)
      y_HD_lim <- trunc((y_HD_check - y_BG_corner) / cellsize)
      mnt_rogne <-
        mnt[(dimension - y_HD_lim):(dimension - y_BG_lim), (x_BG_lim):(x_HD_lim)]
      selection_ligne <- seq(1, dim(mnt_rogne)[1], by = division)
      selection_colonne <- seq(1, dim(mnt_rogne)[2], by = division)
      mnt_rogne <- mnt_rogne[selection_ligne, selection_colonne]
      
      # We create the phi and lambda matrix
      nb_ligne_mnt <- dim(mnt_rogne)[1]
      nb_colonne_mnt <- dim(mnt_rogne)[2]
      M_mnt <- mat.or.vec(nb_ligne_mnt, nb_colonne_mnt)
      phi_mnt <- mat.or.vec(nb_ligne_mnt, nb_colonne_mnt)
      seq1 <- seq(from = y_HD_check,
                  by = -cellsize * division,
                  to = y_BG_check)
      seq1 <- seq1[1:nb_ligne_mnt]
      seq2 <- seq(from = x_BG_check,
                  by = cellsize * division,
                  to = x_HD_check)
      seq2 <- seq2[1:nb_colonne_mnt]
      for (u in 1:nb_colonne_mnt) {
        phi_mnt[, u] <- seq1
      }
      for (u in 1:nb_ligne_mnt) {
        M_mnt[u, ] <- seq2
      }
      
      # We go from EGM altitudes to WGS84 ellipsoid heights
      source('fenetre_info_alti_he.r')
      mnt_rogne <- alti_2_he(seq1, seq2, mnt_rogne, 'EGM96')
      
      # Transformation of matrix into vectors
      M_mnt <- c(t(M_mnt))
      phi_mnt <- c(t(phi_mnt))
      mnt_rogne <- c(t(mnt_rogne))
      
      #
      # Projection : passing ECEF --> ENU (centered on the point projection on the WGS84 ellipsoid)
      #
      pos_R_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_R_xyz)
      pos_R_xy <- c(pos_R_geo[1], pos_R_geo[2])#c(0,0)
      pos_R_geo_approx <-
        c(pos_R_geo[1], atan(R_p_WGS84 / R_0_WGS84 * tan(pos_R_geo[2] * pi / 180)) *
            180 / pi)
      R_local <-
        R_0_WGS84 * (1 - e_WGS84 ^ 2) / sqrt(1 - e_WGS84 ^ 2 * sin(pos_R_geo[2] *
                                                                     pi / 180) * sin(pos_R_geo[2] * pi / 180))
      pas <- cellsize * division  * pi / 180 * R_local * 1.1 # approximation
      f_WGS84 <- R_0_WGS84 / e_WGS84 # flatness of the WGS84 ellipsoid
      lon <-
        M_mnt * pi / 180 # longitude of the considered point in the WGS84 ellipsoid
      lat <-
        phi_mnt * pi / 180 # latitude of the considered point in the WGS84 ellipsoid
      
      N <-
        R_0_WGS84 / (sqrt(1 - e_WGS84 ^ 2 * sin(lat) * sin(lat))) # big normal of the considered point in the WGS84 ellipsoid
      n <- e_WGS84 ^ 2 * N
      alt <- N + mnt_rogne
      long_mnt <- lon * 180 / pi
      lat_mnt <- lat * 180 / pi
      he_mnt <- mnt_rogne
      # Calculation of the WGS84 cartesian coordinates (ECEF)
      X_mnt <- alt * cos(lon) * cos(lat)
      Y_mnt <- alt * sin(lon) * cos(lat)
      Z_mnt <- sin(lat) * (alt - n)
      #  for(u in 1:length(X_mnt)){
      #      result <- pos2geo(R_0_sol,R_p_sol,c(X_mnt,Y_mnt,Z_mnt))
      #      long_mnt[u] <- result[1]
      #      lat_mnt[u] <- result[2]
      #      he_mnt[u] <- result[3]
      #  }
      tkdestroy(ttinfo)
    } # End of the DEM preparation
    
    progressbar <-
      tkProgressBar(
        "Determination of the specular reflection points positions",
        "Progression in %",
        0,
        100,
        0
      ) # Display of a progress bar indicating the calculation progression
    
    
    if (static_receiver) {
      #
      # Initialization of the plot
      #
      
      x11()
      plot(
        x,
        y,
        xlab = 'X_WGS84 (m)',
        ylab = 'Y_WGS84 (m)',
        main = 'Positions of the specular reflection points',
        pch = 11,
        lwd = 2,
        lty = 'dotdash',
        type = 'h',
        col = 'blue',
        xaxs = 'r',
        yaxs = 'r',
        xlim = c(x_min, x_max),
        ylim = c(y_min, y_max)
      )
      
      
      #
      # Creation of a North arrow
      #
      
      pos_Nord <- c(0, 0, R_p_WGS84)
      if ((pos_Rsol_xyz != pos_Nord &
           pos_Rsol_xyz != -pos_Nord)[1]) {
        # If the reception point is NOT the North or South pole
        
        a_droite_SN <-
          (pos_Nord[2] - pos_Rsol_xyz[2]) / (pos_Nord[1] - pos_Rsol_xyz[1])
        b_droite_SN <- pos_Nord[2] - a_droite_SN * pos_Nord[1]
        
        x_fleche <-
          x_centre + sign(pos_Nord[1] - x_centre) * (x_max - x_centre) * 3 / 4 # Abscissa of the final point of the North arrow, which will go from the center of the plot, until the 3/4 maximum.
        y_fleche <-
          a_droite_SN * (x_centre + sign(pos_Nord[1] - x_centre) * (x_max - x_centre) *
                           3 / 4) + b_droite_SN # Ordinate of the final point of the arrow
        arrows(
          x_centre,
          y_centre,
          x_fleche,
          y_fleche,
          angle = 45,
          code = 2,
          lwd = 3
        )  # Creation of the arrow
        text(x_fleche - (x_max - x_min) / 15, y_fleche, 'O', cex = 2) # We add a "N" label at the end of the North arrow
      }
    }
    
    reflexion <- FALSE # boolean who says if the were reflection or not
    sat_visibles <- list()
    increment_sat <- 1
    elev_sat <- list()
    increment_elev <- 1
    
    while (num_ligne < nb_ligne) {
      # We proceed from line to line going through mon_data data, given in input
      num_sat <- mon_data[num_ligne, 1]
      
      #
      # Update of the progress bar
      #
      
      avance <- floor(num_ligne / nb_ligne * 100)
      info <- sprintf("%d%% done", round(avance), 12)
      setTkProgressBar(progressbar,
                       avance,
                       sprintf("Calculation progression (%s)", info),
                       info)
      tcl("update")
      
      
      #
      # We skip a line in the output file between each satellite,  for clarity
      #
      
      ecrit <- c('\n')
      write.table(
        ecrit,
        file = gsub(' ', '', paste(path_output, "output.txt")),
        append = TRUE,
        quote = FALSE,
        sep = " ",
        na = "-",
        row.names = FALSE,
        col.names = FALSE
      )
      if (affichage_output) {
        write.table(
          ecrit,
          file = stdout(),
          append = TRUE,
          quote = FALSE,
          sep = " ",
          na = "-",
          row.names = FALSE,
          col.names = FALSE
        )
      }
      #
      # The current line may either indicate a line, or a satellite (cf .sp3 file pattern)
      #
      
      if (mon_data[num_ligne, 1] == '*') {
        # If the current line indicates a date
        compteur <- compteur + 1
        annee <- as.character(mon_data[num_ligne, 2])
        mois <- as.character(mon_data[num_ligne, 3])
        jour <- as.character(mon_data[num_ligne, 4])
        heure <- as.character(mon_data[num_ligne, 5])
        minute <- as.character(mon_data[num_ligne, 6])
        seconde <-
          as.character(round(as.numeric(mon_data[num_ligne, 7])))
        
        if (nchar(mois) == 1) {
          mois <- gsub(' ', '', paste('0', as.character(mois)))
        }
        if (nchar(jour) == 1) {
          jour <- gsub(' ', '', paste('0', as.character(jour)))
        }
        if (nchar(heure) == 1) {
          heure <- gsub(' ', '', paste('0', as.character(heure)))
        }
        if (nchar(minute) == 1) {
          minute <- gsub(' ', '', paste('0', as.character(minute)))
        }
        if (nchar(seconde) == 1) {
          seconde <- gsub(' ', '', paste('0', as.character(seconde)))
        }
        date_courante <-
          gsub(
            ' ',
            '',
            paste(
              annee,
              '-',
              mois,
              '-',
              jour,
              'T',
              heure,
              ':',
              minute,
              ':',
              seconde,
              'Z'
            )
          )
        
        
        date_jul_cnes <-
          round((as.double(as.POSIXct(
            paste(
              annee,
              '-',
              mois,
              '-',
              jour,
              ' ',
              heure,
              ':',
              minute,
              ':',
              seconde,
              sep = ''
            )
          )) -  as.double(as.POSIXct(
            '1950-01-01 00:00:00'
          ))) / 60 / 60 / 24, 4)
        
        #
        # We write the date in the output file
        #
        
        ecrit <-
          t(
            c(
              ' --------------------------------------------------------------------\n',
              'Date :',
              jour,
              mois,
              annee,
              'at',
              heure,
              'h',
              minute,
              'min et',
              seconde,
              's\n',
              '--------------------------------------------------------------------'
            )
          )
        write.table(
          ecrit,
          file = gsub(' ', '', paste(path_output, "output.txt")),
          append = TRUE,
          quote = FALSE,
          sep = " ",
          na = "-",
          row.names = FALSE,
          col.names = FALSE
        )
        if (affichage_output) {
          write.table(
            ecrit,
            file = stdout(),
            append = TRUE,
            quote = FALSE,
            sep = " ",
            na = "-",
            row.names = FALSE,
            col.names = FALSE
          )
        }
        num_ligne <- num_ligne + 1 # We go to the next line
        
        
        continuer_hauteur <- TRUE
        if (!static_receiver) {
          pos_R_geo <- pos_R_geo_tot[compteur, ]
          pos_R_xyz <- geo2pos(R_0_WGS84, R_p_WGS84, pos_R_geo)
          diff_h <- pos_R_geo[3] - pos_R_geo_tot[1, 3]
          hauteur_recepteur <- hauteur_recepteur_init + diff_h
          if (hauteur_recepteur <= 0) {
            print('Il y a une erreur sur laltitude du recepteur (sense etre sous terre)',
                  quote = FALSE)
            continuer_hauteur <- FALSE
          }
          write.table(
            t(
              c(
                round(pos_R_geo[1], 7),
                round(pos_R_geo[2], 7),
                round(pos_R_geo[3], 3),
                param[12, 2],
                param[13, 2],
                param[14, 2],
                param[15, 2],
                date_courante
              )
            ),
            gsub(' ', '', paste(path_output, "receiver.txt")),
            append = TRUE,
            quote = FALSE,
            sep = " ",
            na = "-",
            row.names = FALSE,
            col.names = FALSE
          )
        }
      }
      
      
      
      else{
        # if the current line indicates a satellite
        
        # 2D (ENU) coordinates of the ground projection of the reception point
        
        pos_E_xyz <-
          c(
            as.numeric(mon_data[num_ligne, 2]) * 1000,
            as.numeric(mon_data[num_ligne, 3]) * 1000,
            as.numeric(mon_data[num_ligne, 4]) * 1000
          ) # cartesian coordinates of the transmittor (satellite) in the IGS system, approximated to the WGS84 system (very small differences)
        pos_E_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_E_xyz)
        
        # Current date
        num_sat <- mon_data[num_ligne, 1]
        
        
        
        # Writing of a file containing the satellites ground paths ('trace_sat_PG01.tmp',trace_sat_PG02.tmp',...)
        nom_fichier_trace_sat <-
          paste("trace_sat_", as.character(num_sat), '.tmp', sep = "")
        write.table(
          t(
            c(
              pos_E_geo[1],
              pos_E_geo[2],
              pos_E_geo[3],
              date_courante,
              as.character(num_sat)
            )
          ),
          gsub(' ', '', paste(
            path_output, nom_fichier_trace_sat
          )),
          append = TRUE,
          quote = FALSE,
          sep = " ",
          na = "-",
          row.names = FALSE,
          col.names = FALSE
        )
        
        if (!static_receiver) {
          #
          # Calculation of the ellipsoid similar to the WGS84 ellipsoid, going through the ground projection of the receiver (= receiver with a zero altitude)
          #
          
          ellipsoide_sol <-
            equ_ellipsoide(R_0_WGS84, R_p_WGS84, pos_R_xyz, hauteur_recepteur) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the ground projection of the receiver
          R_0_sol <- ellipsoide_sol[1]
          R_p_sol <- ellipsoide_sol[2]
          #
          # Calculation of the ellipsoid similar to the WGS84, going through the receiver
          #
          
          ellipsoide_recepteur <-
            equ_ellipsoide(R_0_WGS84, R_p_WGS84, pos_R_xyz, 0) # We get R_0 et R_p (semi-minor axis and semi-major axis) of the ellipsoid "parallel" to the WGS84 ellipsoid, going through the receiver
          R_0_recepteur <- ellipsoide_recepteur[1]
          R_p_recepteur <- ellipsoide_recepteur[2]
          
          #
          # Calculation of the WGS84 cartesian coordinates of the ground projection of the receiver (= receiver with a zero altitude)
          #
          
          pos_Rsol_xyz <- proj_surf(R_0_sol, R_p_sol, pos_R_xyz)
          pos_Rsol_geo <-
            proj_surf(R_0_WGS84, R_p_WGS84, pos_Rsol_xyz)
          pos_Rellips_xyz <-
            proj_surf(R_0_WGS84, R_p_WGS84, pos_R_xyz)
          pos_Rellips_geo <-
            pos2geo(R_0_WGS84, R_p_WGS84, pos_Rellips_xyz)
          
          
          equ_plan_local <-
            equ_plan(R_0_WGS84, R_p_WGS84, pos_R_xyz)
          
          x <- pos_Rsol_xyz[1]
          y <- pos_Rsol_xyz[2]
        }
        
        
        
        # We check whether the satellite is visible with respect to the visibility mask and also the DEM if we take it into account.
        result_visibi_sat <-
          sat_visible(
            pos_E_xyz,
            pos_R_xyz,
            R_0_WGS84,
            R_p_WGS84,
            param[12, 2],
            param[13, 2],
            param[14, 2],
            param[15, 2],
            equ_plan_local
          )
        satellite_visible_masque <-
          result_visibi_sat[[1]] # If the satellite is visible w.r.t the parametered masks
        satellite_visible_mnt <- TRUE
        source('intersection_droite_ellipsoide.r')
        source('calcul_normale.r')
        if (satellite_visible_masque &
            continuer_hauteur) {
          # If the satellite is visible w.r.t the visibility masks, we will check whether it is the case w.r.t the DEM
          if (param[2, 2] == '-') {
            satellite_visible_mnt == TRUE
          } # If the user doesn't want to take a DEM into acocunt, we consider that the DEM does not mask the satellite
          else{
            # If we take a DEM into account
            
            # projection of the receiver on the ellipsoid
            pos_1 <- proj_surf(R_0_sol, R_p_sol, pos_R_xyz)
            
            # projection of the satellite on the ellipsoid
            pos_2 <- proj_surf(R_0_sol, R_p_sol, pos_E_xyz)
            
            increment <- pas / 5
            x_droite <- 0
            y_droite <- 0
            z_droite <- 0
            x_mnt <- 0
            y_mnt <- 0
            z_mnt <- 0
            x_geo <- 0
            y_geo <- 0
            z_geo <- 0
            vec <- (pos_2 - pos_1) / vecnorm(pos_2 - pos_1)
            
            contin <- TRUE
            nombre <- 0
            zer <- 1
            # he_mnt <- he_mnt*0+420.068763716736
            normale <- calcul_normale(pos_1, pos_2)
            
            while (contin) {
              x_droite <- vec[1] * nombre + pos_1[1]
              y_droite <- vec[2] * nombre + pos_1[2]
              z_droite <- vec[3] * nombre + pos_1[3]
              
              # projection of the point on the ellipsoid
              pos <-
                intersection_droite_ellipsoide(
                  R_0_sol,
                  R_p_sol,
                  normale[1],
                  normale[2],
                  normale[3],
                  x_droite,
                  y_droite,
                  z_droite
                )
              # projection of the point on the DEM
              pos_geo <- pos2geo(R_0_sol, R_p_sol, pos)
              if (pos_geo[1] > min(long_mnt) &
                  pos_geo[1] < max(long_mnt) &
                  pos_geo[2] > min(lat_mnt) & pos_geo[2] < max(lat_mnt)) {
                x_geo[zer] <- pos_geo[1]
                y_geo[zer] <- pos_geo[2]
                
                result <-
                  pos2geo(R_0_WGS84,
                          R_p_WGS84,
                          geo2pos(R_0_sol, R_p_sol, c(x_geo[zer], y_geo[zer], 0)))
                x_geo[zer] <- result[1]
                y_geo[zer] <- result[2]
                z_geo[zer] <- result[3]
                
                zer <- zer + 1
                nombre <- nombre + increment
                if (nombre > vecnorm(pos_2 - pos_1)) {
                  contin <- FALSE
                }
              }
              else{
                contin <- FALSE
              }
            }
            
            
            z_interp <-
              interp_mnt(x_geo, y_geo, long_mnt, lat_mnt, he_mnt, param[31, 2])
            #z_interp <- z_geo  # si on veut comparer mnt et ellipsoide
            
            for (zer in 1:length(z_interp)) {
              pos <-
                geo2pos(R_0_WGS84, R_p_WGS84, c(x_geo[zer], y_geo[zer], z_interp[zer]))
              x_mnt[zer] <- pos[1]
              y_mnt[zer] <- pos[2]
              z_mnt[zer] <- pos[3]
            }
            
            #We check if the satellite is visible from the receiver
            satellite_visible_mnt <-
              mnt_visible(x_mnt,
                          y_mnt,
                          z_mnt,
                          pos_R_xyz,
                          pos_E_xyz,
                          R_0_sol,
                          R_p_sol)
          } # end "if we take a DEM into account"
          if (satellite_visible_mnt) {
            # If the satellite is visible w.r.t the DEM --> so, if it is visible w.r.t the DEM AND w.r.t the visibility masks
            if (param[2, 2] != '-') {
              # if we take a DEM into account
              pos_S_3D <-
                reflexion_mnt(
                  x_mnt,
                  y_mnt,
                  z_mnt,
                  pos_R_xyz,
                  pos_E_xyz,
                  as.numeric(param[11, 2]),
                  X_mnt,
                  Y_mnt,
                  Z_mnt,
                  long_mnt,
                  lat_mnt,
                  he_mnt,
                  R_0_sol,
                  R_p_sol
                )
              pos_S <- pos_S_3D
              pos_Svrai_xy_1D <- pos_S_3D[[6]]
              
              if (pos_S_3D[[1]] == FALSE) {
                # if the receiver is not visible from the reflection point, or is the satellite is not visible from the reflection point, or if the calculation did not converge
                ecrit <-
                  t(
                    c(
                      ' Any reflection point could be find (due to the relief, or the chosen tolerance might be too restrictive !)'
                    )
                  )
                write.table(
                  ecrit,
                  file = gsub(' ', '', paste(path_output, "output.txt")),
                  append = TRUE,
                  quote = FALSE,
                  sep = " ",
                  na = "-",
                  row.names = FALSE,
                  col.names = FALSE
                )
                if (affichage_output) {
                  write.table(
                    ecrit,
                    file = stdout(),
                    append = TRUE,
                    quote = FALSE,
                    sep = " ",
                    na = "-",
                    row.names = FALSE,
                    col.names = FALSE
                  )
                }
                pos_S <- list(pos_S_3D[[1]])
              }
              
              if (pos_S_3D[[1]] == TRUE) {
                # if there were reflection
                pos_S_xyz <-
                  c(pos_S_3D[[2]], pos_S_3D[[3]], pos_S_3D[[4]])
                pos_S_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
              }
            }
            else{
              # if we do not take the DEM into account
              
              # Coordinates in the "local" ellipsoid
              pos_Rsol_xyz <- proj_surf(R_0_sol, R_p_sol, pos_R_xyz)
              pos_Rsol_geo <- pos2geo(R_0_sol, R_p_sol, pos_Rsol_xyz)
              
              # Parameters of the osculating sphere
              R_gaussian <-
                R_0_sol ^ 2 * R_p_sol / ((R_0_sol * cos(pos_Rsol_geo[2] * pi / 180)) ^
                                           2 + (R_p_sol * sin(pos_Rsol_geo[2] * pi / 180)) ^ 2)
              pos_centre_geo <-
                c(pos_Rsol_geo[1], pos_Rsol_geo[2], -R_gaussian)
              pos_centre_xyz <-
                geo2pos(R_0_sol, R_p_sol, pos_centre_geo)
              
              # We calculate the satellite elevation
              
              vec_normal <- pos_Rsol_xyz - pos_centre_xyz
              # equation of the normal plane going through the receiver
              a <- vec_normal[1]
              b <- vec_normal[2]
              c1 <- vec_normal[3]
              d <-
                -(a * pos_R_xyz[1] + b * pos_R_xyz[2] + c1 * pos_R_xyz[3])
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
              if (param[17, 2] == 'ellipsoid') {
                # If the user wants to use the ellipsoid algorithm
                
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
                                         as.numeric(param[19, 2]))
                pos_S_xyz <-
                  c(pos_S[[2]][1], pos_S[[3]][1], pos_S[[4]][1])
                pos_S_geo <-
                  pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
              }
              
              else if (param[17, 2] == 'plane') {
                source('reflexion_algo_plan.r')
                pos_R_geo <- pos2geo(R_0_sol, R_p_sol, pos_R_xyz)
                pos_Rfictif_geo <-
                  c(pos_R_geo[1],
                    pos_R_geo[2],
                    pos_R_geo[3] - 2 * hauteur_recepteur)
                pos_Rfictif_xyz <-
                  geo2pos(R_0_sol, R_p_sol, pos_Rfictif_geo)
                pos_Esol_xyz <- proj_plan(c(a, b, c1, d), pos_E_xyz)
                pos_S <-
                  reflexion_algo_plan(
                    pos_E_xyz,
                    pos_Esol_xyz,
                    pos_R_xyz,
                    pos_Rsol_xyz,
                    pos_Rfictif_xyz,
                    hauteur_recepteur,
                    distance_plan,
                    a,
                    b,
                    c1,
                    d
                  ) # We calculate the reflection point position
                pos_S_xyz <-
                  c(pos_S[[2]][1], pos_S[[3]][1], pos_S[[4]][1])
                pos_S_geo <-
                  pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
              }
              else if (param[17, 2] == 'sphere') {
                # If the user wants to use the sphere algorithm
                
                pos_S <-
                  reflexion_algo_helm(
                    R_gaussian,
                    pos_centre_xyz,
                    pos_E_xyz,
                    pos_R_xyz,
                    elevation_satellite,
                    hauteur_recepteur,
                    as.numeric(param[20, 2])
                  )
                pos_S_xyz <-
                  c(pos_S[[2]][1], pos_S[[3]][1], pos_S[[4]][1])
                pos_S_geo <-
                  pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
                R_sphere <- as.numeric(pos_S[7])
              }
              
              
            }
            
            #
            ######### Tropospheric correction #################
            #
            grub <- 1
            ab <- 0
            bb <- 0
            cb <- 0
            
            if (correction_tropo) {
              # if we want to take the tropo correction into account
              
              # we load the necessary functions
              
              source('calc_amfx.r')
              source('amfx_prlx.r')
              
              # We take the name of the AMF file, calculated until the end of the troposphere (6450 km)
              nom_fichier_amf <-
                paste(path_input_amfx, nom_fichier_amfx_sat, sep = '') #paste("AMF-X-F3A4C5G2-COR1_COR2_COR3_LEM1_LEM2-22921.3750_23130.5000.amf",sep='/')
              # We take the name of the station we want to get information in the AMF file (COR1, LEM1, etc...)
              station <- nom_station_amfx
              
              ellipsoide_S <-
                equ_ellipsoide(R_0_WGS84, R_p_WGS84, pos_S_xyz, 0) # We get R_0 and R_p (semi-major and -minor axis of the ellipsoid)
              R_0_S <- ellipsoide_S[1]
              R_p_S <- ellipsoide_S[2]
              equ_plan_local_S <-
                equ_plan(R_0_WGS84, R_p_WGS84, pos_S_xyz) # R_0_S,R_p_S
              result <-
                sat_visible(
                  pos_E_xyz,
                  pos_S_xyz,
                  R_0_WGS84,
                  R_p_WGS84,
                  0,
                  90,
                  0,
                  360,
                  equ_plan_local_S
                )
              pos_R_geo <-
                pos2geo(R_0_WGS84, R_p_WGS84, pos_R_xyz)
              # End of the calculation of the satellite elevation and azimuth, seen from the reflection point
              pos_S_geo <-
                pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
              E_sat <-
                result[[3]] # satellite elevation seen from the reflection point (rad)
              pos_s <- pos_S_xyz
              pos_S_projete <-   pos_S_xyz
              if (param[2, 2] != '-') {
                pos_S_xy_1D <- pos_Svrai_xy_1D
              }
              #x0 <- pos_R_xy_1D[1]
              #                              x1 <- pos_Stemp_xy[1]
              # liste_num <- #seq(from=pos_Stemp_xy[1],to=pos_R_xy_1D[1],by=-as.numeric(param[31,2]))
              for (nb_it in 1:2) {
                # one iteration would almost be sufficient ...
                
                
                #
                # pos_S_xyz is a vector containing the 3 cartesian coordinates of the reflection point, calculated without taking the tropo correction into account
                # pos_E_xyz is a vector containing the 3 cartesian coordinates of the satellite
                # pos_R_xyz is a vector containing the 3 cartesian coordinates of the receiver
                # pos_S_geo, pos_E_geo et pos_R_geo are vectors containing the 3 geodesic coordinates
                
                #
                # 1. We first consider the part of the reflected signal going from the satellite to the reflection point --> we want to determine E_sit, the real elevation of the wave coming to the reflection point
                #
                
                # We first calculate the satellite elevation and azimuth, seen from that reflection point (i.e. without taking the tropo into account)
                
                ellipsoide_S <-
                  equ_ellipsoide(R_0_WGS84, R_p_WGS84, pos_S_xyz, 0) # We get R_0 et R_p (semi-major axis and semi-minor axis of the WGS84 ellipsoid)
                R_0_S <- ellipsoide_S[1]
                R_p_S <- ellipsoide_S[2]
                equ_plan_local_S <-
                  equ_plan(R_0_WGS84, R_p_WGS84, pos_S_xyz) # R_0_S,R_p_S
                result <-
                  sat_visible(
                    pos_E_xyz,
                    pos_S_xyz,
                    R_0_WGS84,
                    R_p_WGS84,
                    0,
                    90,
                    0,
                    360,
                    equ_plan_local_S
                  )
                pos_R_geo <-
                  pos2geo(R_0_WGS84, R_p_WGS84, pos_R_xyz)
                
                
                E_sat <-
                  result[[3]] # satellite elevation seen from the reflection point (rad)
                A_sat <-
                  (pi + result[[2]]) %% (2 * pi)# satellite azimuth seen from the reflection point (rad)  !!!! Beware --> we do 180o + azimuth to have the azimuth with respect to the South (which is mandatory to compute the AMF --> cf Pascal GEGOUT article)
                R_sat <-
                  sqrt(pos_E_xyz[1] ^ 2 + pos_E_xyz[2] ^ 2 + pos_E_xyz[3] ^ 2)  # satellite radius (geocentric coordinates)
                T_sat <-
                  pi / 2 - (asin(pos_E_xyz[3] / R_sat)) # geocentric colatitude of the satellite(rad)
                L_sat <-
                  pos_E_geo[1] * pi / 180# satellite longitude (rad)
                R_sit <-
                  vecnorm(pos_R_xyz) - hauteur_recepteur#sqrt(pos_S_xyz[1]^2+pos_S_xyz[2]^2+pos_S_xyz[3]^2)# reflection point radius (geocentric coordinates)
                T_sit <-
                  pi / 2 - (asin(pos_S_xyz[3] / R_sit))# geocentric colatitude of the reflection point (rad)
                F_sit <-
                  pos_S_geo[2] * pi / 180# "geodetic" latitude of the reflection point(rad)
                L_sit <-
                  pos_S_geo[1] * pi / 180 # longitude of the reflection point (rad)
                R_sta <-
                  sqrt(pos_R_xyz[1] ^ 2 + pos_R_xyz[2] ^ 2 + pos_R_xyz[3] ^ 2) # receiver radius (geocentric coordinates)
                R_top <- 6450000 # tropo radius
                T_sta <-
                  pi / 2 - (asin(pos_R_xyz[3] / R_sta))
                L_sta <- pos_R_geo[1] * pi / 180
                
                # We solve the parallax to get E_sit, the real elevation to which comes the wave on the reflection point
                
                calcu_para <-
                  calc_prlx(
                    nom_fichier_amf,
                    station,
                    date_jul_cnes,
                    E_sat,
                    A_sat,
                    T_sat,
                    L_sat,
                    R_sat,
                    T_sit,
                    F_sit,
                    R_sit,
                    L_sit,
                    R_top
                  )
                E_sit <- calcu_para[1]
                G_top <- calcu_para[3]
                A_sat <- calcu_para[5]
                eee <- calcu_para[8]
                D_ray <- calcu_para[7]
                
                #
                #
                # 2. We now consider the part of the reflected signal going from the reflection point to the receiver
                #
                #
                
                #  We chant the AMF file (we take the one calculated until the receiver and not the one calculated until the end of the tropo)
                nom_fichier_amf <-
                  paste(path_input_amfx, nom_fichier_amfx_sta, sep = '') #paste("AMF-X-F3A4C5G2-COR4_LEM4-22921.3750_22923.2500.amf",sep='/')
                
                A_sta <-
                  (A_sat + pi) %% (2 * pi) # New azimuth ==> equal to the previous one + 180 o
                R_top <-
                  R_sta # Now, we consider that the tropo ends at the receiver height
                # We calculate G_top, the staight distance (== geometric) between the reflection point and the receiver
                E_cor <- (pi / 2 - T_sit - F_sit) * cos(A_sta)
                A_cor <-
                  (pi / 2 - T_sit - F_sit) * sin(A_sta) * tan(E_sit)
                D_top <-
                  E_sit - calc_amfx(nom_fichier_amf,
                                    1,
                                    station,
                                    date_jul_cnes,
                                    E_sit + E_cor,
                                    A_sta + A_cor) * pi / 180
                P_top <-
                  acos(R_sit / R_top * cos(D_top)) - D_top
                G_top <-
                  sqrt(R_sit ^ 2 + R_top ^ 2 - 2 * R_sit * R_top * cos(P_top))
                
                T_sit <-
                  pi / 2 - (asin(pos_S_xyz[3] / R_sit))
                L_sit <- pos_S_geo[1] * pi / 180
                CT_sit <- cos(T_sit)
                ST_sit <- sin(T_sit)
                CT_sta <- cos(T_sta)
                ST_sta <- sin(T_sta)
                CP_sta <-
                  CT_sit * CT_sta + ST_sit * ST_sta * cos(L_sit - L_sta)
                P_sta <- acos(CP_sta)
                SP_sta <- sin(P_sta)
                D_sta <-
                  R_sta * sqrt(1 + (R_sit / R_sta) * (R_sit / R_sta - 2 * CP_sta))
                CE_sta <- R_sta / D_sta * SP_sta
                SE_sta <-
                  (R_sta ^ 2 - R_sit ^ 2 - D_sta ^ 2) / (2 * R_sit * D_sta)
                E_sta <- sign(SE_sta) * acos(CE_sta)
                pos_S_projete <- pos_S_xyz
                
                ecart <- E_sit - eee
                
                #if(ecart<0){sdfsefse}
                #                              while(ecart > 0.00000001){
                #                                facteur <- 1/100
                #                                vec_S_R <- (pos_R_xyz-pos_S_projete)/vecnorm(pos_R_xyz-pos_S_projete)
                #                                pos_S_non_projete <- pos_S_projete + facteur*vec_S_R
                #                                pos_S_projete <- proj_surf(R_sphere,R_sphere,pos_S_non_projete)
                #                                pos_S_geo_projete <- pos2geo(R_0_WGS84,R_p_WGS84,pos_S_projete)
                #                                T_sit <- pi/2-(asin(pos_S_projete[3]/R_sit))
                #                                L_sit <- pos_S_geo_projete[1]*pi/180
                #                                CT_sit <- cos(T_sit)
                #                                ST_sit <- sin(T_sit)
                #                                CT_sta <- cos(T_sta)
                #                                ST_sta <- sin(T_sta)
                #                                CP_sta <- CT_sit*CT_sta+ST_sit*ST_sta*cos(L_sit-L_sta)
                #                                 P_sta <- acos(CP_sta)
                #                                 SP_sta <- sin(P_sta)
                #                                D_sta <- R_sta*sqrt(1+(R_sit/R_sta)*(R_sit/R_sta-2*CP_sta))
                #                                CE_sta <- R_sta/D_sta*SP_sta
                #                                SE_sta <- (R_sta^2-R_sit^2-D_sta^2)/(2*R_sit*D_sta)
                #                                E_sta <- sign(SE_sta)*acos(CE_sta)
                #
                #                              ecart <-  E_sit - E_sta
                #
                #
                #                              }
                #
                
                if (param[2, 2] == '-') {
                  # if we do not integrate a DEM
                  
                  iteration <- 1
                  ecart <- 12
                  while (ecart >= 0.0005) {
                    # While the angle we find is greater than the fixed tolerance
                    if (ecart > 1) {
                      facteur <- 1 / 10
                    }
                    else{
                      facteur <- 1 / 100
                    }
                    vec_S_R <-
                      (pos_R_xyz - pos_S_projete) / vecnorm(pos_R_xyz - pos_S_projete)
                    pos_S_non_projete <-
                      pos_S_projete + facteur * vec_S_R
                    if (param[17, 2] == 'sphere') {
                      pos_S_projete <- proj_surf(R_sphere, R_sphere, pos_S_non_projete)
                    }
                    else if (param[17, 2] == 'ellipsoid') {
                      pos_S_projete <- proj_surf(R_0_sol, R_p_sol, pos_S_non_projete)
                    }
                    ecart <-
                      vecnorm(pos_R_xyz - pos_S_projete) - G_top
                    iteration <- iteration + 1
                    
                    if (iteration > 10000) {
                      print('Too many iterations !!! Check why !')
                      break
                    }    # If there are too many iterations, we go out of the loop
                  }
                  
                  pos_S_xyz <-
                    pos_S_projete # We get the new position of the specular reflection point
                  ### All that follows in commentary, is if we want to solve the parallax with the new position of the reflection point and to check if E_sit is effectively equal in both sides of the reflection. In practice, one or two iterations are needed to be the case.
                  ### ==> That's why we have to iterate all the previous step, because the new position of the reflection point check the good distance from the receiver (G_top), but does not necessarely check E_sit(sat)=E_sit(rocepteur) any more, because, the E_sit was calculated for the initial position of the reflection point.
                  #E_sat <- result[[3]] # satellite elevation (rad)
                  #              A_sat <- result[[2]]# satellite azimuth(rad)
                  #              R_sat <- sqrt(pos_E_xyz[1]^2+pos_E_xyz[2]^2+pos_E_xyz[3]^2)  # satellite radius (geocentric coordinates)
                  #              T_sat <- pi/2-(asin(pos_E_xyz[3]/R_sat)) # geocentric colatitude of the satellite (rad)
                  #              L_sat <- pos_E_geo[1]*pi/180 # Satellite longitude(rad)
                  #              R_sit <- sqrt(pos_S_xyz[1]^2+pos_S_xyz[2]^2+pos_S_xyz[3]^2)# site radius (geocentric coordinates)
                  #              T_sit <- pi/2-(asin(pos_S_xyz[3]/R_sit))# geocentric colatitude of the site (rad)
                  #              F_sit <- pos_S_geo[2]*pi/180# geodetic latitude of the site (rad)
                  #              L_sit <- pos_S_geo[1]*pi/180 # site longitude (rad)
                  #              R_sta <- sqrt(pos_R_xyz[1]^2+pos_R_xyz[2]^2+pos_R_xyz[3]^2)
                  #              R_top <- 6450000
                  #
                  #              E_sit <- calc_prlx(nom_fichier_amf,station,date_jul_cnes,E_sat,A_sat,T_sat,L_sat,R_sat,T_sit,F_sit,R_sit,L_sit,R_top)[1]
                  #              print(E_sit)
                }
                
                else if (param[2, 2] != '-') {
                  # if we integrate a DEM
                  
                  #  We calculate a fictive position of the DEM
                  #  pos_Stemp_xy_1D
                  #  pos_E_xy_1D
                  #  pos_R_xy_1D
                  pos_Efictif_xy_1D <- c(0, 0)
                  pos_Efictif_xy_1D[1] <-
                    pos_S_xy_1D[1] + vecnorm(pos_S_xy_1D - pos_E_xy_1D) * sin(pi / 2 - E_sit)
                  pos_Efictif_xy_1D[2] <-
                    pos_S_xy_1D[2] + vecnorm(pos_S_xy_1D - pos_E_xy_1D) * cos(pi / 2 - E_sit)
                  # We pass from 1D fictive satellite coordinates to 2D (ENU)
                  resultat <-
                    passage_1D_ENU(
                      pos_Efictif_xy_1D[1],
                      pos_Efictif_xy_1D[2],
                      pos_E_xy_1D,
                      pos_R_xy_1D,
                      pos_R_xyz_mnt,
                      pos_E_xyz_mnt
                    )
                  # We pass from the 2D (ENU) reflection point coordinates to the 3D (ECEF) coordinates
                  pos_Efictif_3D <-
                    ENU_2_ECEF(
                      pos_Rellips_geo,
                      pos_Rellips_xyz,
                      c(resultat[[1]], resultat[[2]], resultat[[3]])
                    )
                  pos_E_xyz <-
                    c(pos_Efictif_3D[[1]],
                      pos_Efictif_3D[[2]],
                      pos_Efictif_3D[[3]])
                  pos_E_xy_1D <- pos_Efictif_xy_1D
                  resultat <-
                    ECEF_2_ENU(
                      pos_Rellips_geo,
                      pos_Rellips_xyz,
                      pos_E_xyz[1],
                      pos_E_xyz[2],
                      pos_E_xyz[3]
                    )
                  pos_E_xyz_mnt <-
                    c(resultat[[1]], resultat[[2]], resultat[[3]])
                  # We transform the 2D (ENU) coordinates to 1D coordinates
                  resultat <-
                    passage_2D_1D(
                      X_mnt,
                      Y_mnt,
                      mnt_rogne,
                      pos_E_xyz_mnt,
                      pos_R_xyz_mnt,
                      pas,
                      hauteur_recepteur,
                      pos_Rsol_z
                    )
                  X_1D <- resultat[[1]]
                  mnt_1D <- resultat[[2]]
                  
                  # We check wether the satellite is visible w.r.t the visibility masks and w.r.t the DEM we take into account
                  result_visibi_sat <-
                    sat_visible(
                      pos_E_xyz,
                      pos_R_xyz,
                      R_0_WGS84,
                      R_p_WGS84,
                      param[12, 2],
                      param[13, 2],
                      param[14, 2],
                      param[15, 2],
                      equ_plan_local
                    )
                  satellite_visible_masque <-
                    result_visibi_sat[[1]] # if the satellite is visible or not w.r.t the visibility masks
                  satellite_visible_mnt <- TRUE
                  if (satellite_visible_masque) {
                    # If the satellite is visibile w.r.t the visibility masks, we will check if it is visible w.r.t the DEM
                    # 2D (ENU) satellite coordinates
                    resultat <-
                      ECEF_2_ENU(
                        pos_Rellips_geo,
                        pos_Rellips_xyz,
                        pos_E_xyz[1],
                        pos_E_xyz[2],
                        pos_E_xyz[3]
                      )
                    pos_E_xyz_mnt <-
                      c(resultat[[1]], resultat[[2]], resultat[[3]])
                    satellite_visible_mnt <-
                      mnt_visible_2D(
                        mnt_rogne,
                        X_mnt,
                        Y_mnt,
                        pos_E_xyz_mnt,
                        pos_R_xyz_mnt,
                        pas,
                        dimension
                      )
                    if (satellite_visible_mnt) {
                      # if the satellite is not visible w.r.t the DEM --> so, if it is not visible w.r.t. the DEM AND w.r.t the visibility masks
                      # We transform the 2D (ENU) coordinates to 1D coordinates
                      # resultat <- passage_2D_1D(X_mnt,Y_mnt,mnt_rogne,pos_E_xyz_mnt,pos_R_xyz_mnt,pas,hauteur_recepteur,pos_Rsol_z)
                      #         # X_1D <- resultat[[1]]
                      ##          mnt_1D <- resultat[[2]]
                      ##
                      ##          pos_E_xy_1D <- resultat[[3]]
                      ##          pos_R_xy_1D <- resultat[[4]]
                      
                      # We calculate the reflection point position on the ellipsoid (sphere algorithm)
                      distance_plan <-
                        dist_plan(equ_plan_local, pos_E_xyz) # We calculate the distance between the satellite and its projection on the plane going through the receiver, tangent to the ellipsoid going through the receiver
                      elevation_satellite <-
                        asin(distance_plan / vecnorm(pos_E_xyz - pos_R_xyz))
                      pos_S_ellipsoide <-
                        reflexion_algo_helm(
                          R_0_sol,
                          R_p_sol,
                          pos_E_xyz,
                          pos_R_xyz,
                          elevation_satellite,
                          hauteur_recepteur,
                          as.numeric(param[20, 2])
                        ) # We calculate the reflection point position
                      pos_S_x_1D <-
                        distance_entre_points(
                          pos_R_xyz,
                          c(
                            pos_S_ellipsoide[[2]],
                            pos_S_ellipsoide[[3]],
                            pos_S_ellipsoide[[4]]
                          )
                        )
                      # We calculate the reflection point position taking the DEM into account
                      pos_S_1D <-
                        reflexion_mnt(
                          mnt_1D,
                          X_1D,
                          pos_E_xy_1D,
                          pos_S_x_1D,
                          pos_R_xy_1D,
                          as.numeric(param[11, 2])
                        )
                      pos_Svrai_xy_1D <- pos_S_1D[[6]]
                      
                      if (pos_S_1D[[1]] == FALSE) {
                        # if the receiver is not visible from the reflection point, or if the satellite is not visible from the reflection point, or if the reflection calculation did not converge
                        ecrit <-
                          t(
                            c(
                              ' Any specular reflection point could be find (due to the relief or choosen tolerance too restrictive)'
                            )
                          )
                        write.table(
                          ecrit,
                          file = gsub(
                            ' ',
                            '',
                            paste(path_output, "output.txt")
                          ),
                          append = TRUE,
                          quote = FALSE,
                          sep = " ",
                          na = "-",
                          row.names = FALSE,
                          col.names = FALSE
                        )
                        if (affichage_output) {
                          write.table(
                            ecrit,
                            file = stdout(),
                            append = TRUE,
                            quote = FALSE,
                            sep = " ",
                            na = "-",
                            row.names = FALSE,
                            col.names = FALSE
                          )
                        }
                        pos_S <- list(pos_S_1D[[1]])
                      }
                      
                      
                      # We transform the reflection point coordinates from 1D to 2D (ENU)
                      resultat <-
                        passage_1D_ENU(
                          pos_S_1D[[2]],
                          pos_S_1D[[3]],
                          pos_E_xy_1D,
                          pos_R_xy_1D,
                          pos_R_xyz_mnt,
                          pos_E_xyz_mnt
                        )
                      
                      pos_S_xy_1D <- pos_Svrai_xy_1D
                      # We transform the reflection point coordinates from 2D (ENU) to 3D (ECEF)
                      pos_S_3D <-
                        ENU_2_ECEF(
                          pos_Rellips_geo,
                          pos_Rellips_xyz,
                          c(resultat[[1]], resultat[[2]], resultat[[3]])
                        )
                      pos_S <-
                        list(
                          pos_S_1D[[1]],
                          pos_S_3D[[1]],
                          pos_S_3D[[2]],
                          pos_S_3D[[3]],
                          pos_S_1D[[4]],
                          pos_S_1D[[5]]
                        )
                      pos_S_xyz <-
                        c(pos_S_3D[[1]], pos_S_3D[[2]], pos_S_3D[[3]])
                      pos_S_geo <-
                        pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
                    }
                  }
                  
                  else{
                    pos_S <- list(FALSE, 1, 1, 1, 1, 1)
                  }
                  
                  #    iteration <- 1
                  #                    ecart <- 100
                  #
                  #                                plot(c(0,0,100000,100000),c(100,0,100,0))
                  #                    while(ecart > 0.05) { # Tant que l'angle que l'on trouve est suporieur o la tolorance fixoe
                  #                            if(iteration > 10000000){break} # S'il y a trop d'itorations, on sort de la boucle
                  #                            x_courant <- x0+iteration/1000
                  #                            if(is.na(x_courant)){break}
                  #                            pos_Stemp_xy <- c(x_courant,interpolation_1D(x_courant,mnt_1D,X_1D))
                  #                            ecart <- vecnorm(pos_R_xy_1D-pos_Stemp_xy)-G_top
                  #                             points(iteration,ecart)
                  #                             print(ecart)
                  #                            if(ecart > 50000){iteration <- iteration + 500}
                  #                            if(x_courant>2*x1){dsfsefsef}
                  #                            #if(ecart > 1){iteration <- iteration + 50}
                  ##                            else if(ecart > 0.3 & ecart < 1){iteration <- iteration + 10}
                  #                            else {iteration <- iteration + 1}
                  #                    }
                  #
                  #                     pos_S_xy <- pos_Stemp_xy
                  
                }
                
                
                
                
                #print(pos_S_xyz)
              } # end of the iterative loop
              
              
              #             if(param[2,2] != '-'){
              #
              #                    # On transforme les coordonnoes du point de roflexion de la 1D vers la 2D (ENU)
              #                    resultat <- passage_1D_ENU(pos_S_xy[1],pos_S_xy[2],pos_E_xy_1D,pos_R_xy_1D,pos_R_xyz_mnt,pos_E_xyz_mnt)
              #
              #                    # On transforme les coordonnoes du point de roflexion de la 2D (ENU) vers la 3D (ECEF)
              #                    pos_S_3D <- ENU_2_ECEF(pos_Rellips_geo,pos_Rellips_xyz,c(resultat[[1]],resultat[[2]],resultat[[3]]))
              #                    pos_S <- list(TRUE,pos_S_3D[[1]],pos_S_3D[[2]],pos_S_3D[[3]],E_sit,0)
              #                    pos_S_xyz <- c(pos_S_3D[[1]],pos_S_3D[[2]],pos_S_3D[[3]])
              #            }
              
              #pos_S_xyz <- c(pos_S[[2]][1],pos_S[[3]][1],pos_S[[4]][1])
              
              #print(pos_S_geo_1)
              
              
              #
            }
            
            # print(vecnorm(azeaz-pos_S_xyz))
            
            
            #
            #
            ##### End of tropo correction
            #
            #
            #pos_S_geo <- pos2geo(R_0_WGS84,R_p_WGS84,pos_S_xyz)
            #        write.table(pos_S_xyz,'tropo_xyz.txt',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
            #         write.table(pos_S_geo,'tropo_geo.txt',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
            #
            if (pos_S[[1]]) {
              # If a reflection point was found
              
              
              reflexion <- TRUE
              for (num_pt_reflexion in 1:1) {
                #length(pos_S[[2]])){
                distance_plan <-
                  dist_plan(equ_plan_local, pos_E_xyz) # We calculate the distance between the satellite and its projection on the plane going through the receiver, tangent to the ellipsoid going through the receiver
                elevation_satellite <-
                  asin(distance_plan / vecnorm(pos_E_xyz - pos_R_xyz))
                
                if (param[2, 2] != '-') {
                  pos_S_xyz <-
                    c(pos_S[[2]][num_pt_reflexion], pos_S[[3]][num_pt_reflexion], pos_S[[4]][num_pt_reflexion])
                  pos_S_geo <-  pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
                }
                x[num_ligne] <- pos_S_xyz[1]
                y[num_ligne] <- pos_S_xyz[2]
                
                pos_S_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_S_xyz)
                M[iter] <-
                  pos_S_geo[1] * pi / 180 # Longitude of the reflection point
                L[iter] <-
                  pos_S_geo[2] * pi / 180 # Latitude of the reflection point
                
                #
                # We write the results in the output.txt and pts_reflexion.txt files
                #
                
                ecrit <-
                  t(
                    c(
                      ' The satellite',
                      as.character(mon_data[num_ligne, 1]),
                      'is not visible from the receiver.\n',
                      'The calculation of the specular reflection point converged. The reflection angle is:',
                      pos_S[[5]][num_pt_reflexion],
                      'degrees. (difference between the incidence angle and reception angle (depends on the convergence):',
                      pos_S[[6]][num_pt_reflexion],
                      'degrees).'
                    )
                  )
                write.table(
                  ecrit,
                  file = gsub(' ', '', paste(path_output, "output.txt")),
                  append = TRUE,
                  quote = FALSE,
                  sep = " ",
                  na = "-",
                  row.names = FALSE,
                  col.names = FALSE
                )
                if (affichage_output) {
                  write.table(
                    ecrit,
                    file = stdout(),
                    append = TRUE,
                    quote = FALSE,
                    sep = " ",
                    na = "-",
                    row.names = FALSE,
                    col.names = FALSE
                  )
                }
                
                # Fresnel surface
                result_fresnel <-
                  surf_Fresnel(longueur_onde,
                               hauteur_recepteur,
                               elevation_satellite)
                
                if (dessiner_fresnel == TRUE) {
                  # May S(Xo,Yo,Zo) be the first Fresnel ellipse center, A(Xa,Ya,Za) be the farther point of the ellipse from the center (and AA its opposite), B(Xb,Yb,Zb) be the closest point of the ellipse from the center (and BB its opposite)
                  equ_plan_local_S <-
                    equ_plan(R_0_WGS84, R_p_WGS84, pos_S_xyz)
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
                    equ_ellipsoide(R_0_WGS84,
                                   R_p_WGS84,
                                   coord_AA,
                                   result_fresnel[2]) # We get R_0 and R_p (semi-major and -minor axis of the ellipsoid parallel to the WGS84 and going through the ground point
                  R_0_b <- ellipsoide_b[1]
                  R_p_b <- ellipsoide_b[2]
                  pos_AA_BB_xyz <- proj_surf(R_0_b, R_p_b, coord_AA)
                  coord_AA_BB <- pos_AA_BB_xyz
                  
                }# end if(dessiner_fresnel)
                
                # Current date
                date_courante <-
                  gsub(
                    ' ',
                    '',
                    paste(
                      annee,
                      '-',
                      mois,
                      '-',
                      jour,
                      'T',
                      heure,
                      ':',
                      minute,
                      ':',
                      seconde,
                      'Z'
                    )
                  )
                
                # Creation of a text file containing all the specular reflection points positions  --> 'specular_points_galileo.txt' or 'specular_points_glonass.txt' or 'specular_points_gps.txt'
                if (substr(as.character(num_sat), 1, 2) == 'PR')
                {
                  nom_fichier_pts_reflexion <- nom_fichier_pts_reflexion_glonass
                  nom_fichier_sat <- nom_fichier_sat_glonass
                }
                else if (substr(as.character(num_sat), 1, 2) == 'PG')
                {
                  nom_fichier_pts_reflexion <- nom_fichier_pts_reflexion_gps
                  nom_fichier_sat <- nom_fichier_sat_gps
                }
                
                else if (substr(as.character(num_sat), 1, 2) == 'PE')
                {
                  nom_fichier_pts_reflexion <- nom_fichier_pts_reflexion_galileo
                  nom_fichier_sat <- nom_fichier_pts_reflexion_galileo
                }
                if (dessiner_fresnel) {
                  write.table(
                    t(
                      c(
                        round(pos_S_geo[1], 7),
                        round(pos_S_geo[2], 7),
                        round(pos_S_geo[3], 3),
                        date_courante,
                        iter,
                        as.character(num_sat),
                        result_visibi_sat[[2]] * 180 / pi,
                        result_visibi_sat[[3]] * 180 / pi,
                        result_fresnel,
                        coord_A,
                        coord_B,
                        coord_AA_BB
                      )
                    ),
                    gsub(
                      ' ',
                      '',
                      paste(path_output, nom_fichier_pts_reflexion)
                    ),
                    append = TRUE,
                    quote = FALSE,
                    sep = " ",
                    na = "-",
                    row.names = FALSE,
                    col.names = FALSE
                  )
                }
                else{
                  write.table(
                    t(
                      c(
                        round(pos_S_geo[1], 7),
                        round(pos_S_geo[2], 7),
                        round(pos_S_geo[3], 3),
                        date_courante,
                        iter,
                        as.character(num_sat),
                        result_visibi_sat[[2]] * 180 / pi,
                        result_visibi_sat[[3]] * 180 / pi,
                        result_fresnel
                      )
                    ),
                    gsub(
                      ' ',
                      '',
                      paste(path_output, nom_fichier_pts_reflexion)
                    ),
                    append = TRUE,
                    quote = FALSE,
                    sep = " ",
                    na = "-",
                    row.names = FALSE,
                    col.names = FALSE
                  )
                }
                
                # Creation of a text file containing all the satellites positions  --> 'satellites_gps.txt' or 'satellites_glonass.txt' or 'satellites_galileo.txt'
                pos_E_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_E_xyz)
                write.table(
                  t(
                    c(
                      round(pos_E_geo[1], 7),
                      round(pos_E_geo[2], 7),
                      round(pos_E_geo[3], 3),
                      date_courante,
                      as.character(num_sat),
                      result_visibi_sat[[2]] * 180 / pi,
                      result_visibi_sat[[3]] * 180 / pi,
                      iter
                    )
                  ),
                  gsub(' ', '', paste(
                    path_output, nom_fichier_sat
                  )),
                  append = TRUE,
                  quote = FALSE,
                  sep = " ",
                  na = "-",
                  row.names = FALSE,
                  col.names = FALSE
                )
                
                # Creation of a temporary file which will be used to display the reflected waves trajectories in the KML files --> 'traj_ondes_reflechies.tmp'
                pos_R_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_R_xyz)
                matrice_ecrit <- mat.or.vec(3, 6)
                matrice_ecrit[1, ] <-
                  c(pos_E_geo,
                    date_courante,
                    iter,
                    as.character(num_sat))
                matrice_ecrit[2, ] <-
                  c(pos_S_geo,
                    date_courante,
                    iter,
                    as.character(num_sat))
                matrice_ecrit[3, ] <-
                  c(pos_R_geo,
                    date_courante,
                    iter,
                    as.character(num_sat))
                nom_fichier_temp_traj_ondes_reflechies <-
                  'traj_ondes_reflechies.tmp'
                write.table(
                  matrice_ecrit,
                  gsub(
                    ' ',
                    '',
                    paste(
                      path_output,
                      nom_fichier_temp_traj_ondes_reflechies
                    )
                  ),
                  append = TRUE,
                  quote = FALSE,
                  sep = " ",
                  na = "-",
                  row.names = FALSE,
                  col.names = FALSE
                )
                
                # Creation of a temporary file which will be used to display the direct waves trajectories in the KML files --> 'traj_ondes_direct.tmp'
                matrice_ecrit <- mat.or.vec(2, 6)
                matrice_ecrit[1, ] <-
                  c(pos_E_geo,
                    date_courante,
                    iter,
                    as.character(num_sat))
                matrice_ecrit[2, ] <-
                  c(pos_R_geo,
                    date_courante,
                    iter,
                    as.character(num_sat))
                nom_fichier_temp_traj_ondes_directes <-
                  'traj_ondes_directes.tmp'
                write.table(
                  matrice_ecrit,
                  gsub(
                    ' ',
                    '',
                    paste(
                      path_output,
                      nom_fichier_temp_traj_ondes_directes
                    )
                  ),
                  append = TRUE,
                  quote = FALSE,
                  sep = " ",
                  na = "-",
                  row.names = FALSE,
                  col.names = FALSE
                )
                iter <-
                  iter + 1 # Number of visible satellites whose calculations converged
                
                if (static_receiver) {
                  #
                  # Plot of the points in the graph
                  #
                  points(
                    x[num_ligne],
                    y[num_ligne],
                    col = (num_ligne - compteur) %% nombre_satellites + 1,
                    pch = 16
                  )  # The color is defined in this way, in order to that the reflection points of a same satellite are in the same color. There is +1 to avoid to put points in white.
                }
              }
            }
            else{
              # If the calculation of the specular reflection point did not converge
              ecrit <-
                t(
                  c(
                    ' The reflection of the wave emitted by the satellite ',
                    as.character(mon_data[num_ligne, 1]),
                    ' is not possible'
                  )
                )
              write.table(
                ecrit,
                file = gsub(' ', '', paste(path_output, "output.txt")),
                append = TRUE,
                quote = FALSE,
                sep = " ",
                na = "-",
                row.names = FALSE,
                col.names = FALSE
              )
              if (affichage_output) {
                write.table(
                  ecrit,
                  file = stdout(),
                  append = TRUE,
                  quote = FALSE,
                  sep = " ",
                  na = "-",
                  row.names = FALSE,
                  col.names = FALSE
                )
              }
            }
          }
          else{
            # If the DEM hides the satellite
            ecrit <-
              t(c(
                ' The relief hides the satellite',
                as.character(mon_data[num_ligne, 1])
              ))
            write.table(
              ecrit,
              file = gsub(' ', '', paste(path_output, "output.txt")),
              append = TRUE,
              quote = FALSE,
              sep = " ",
              na = "-",
              row.names = FALSE,
              col.names = FALSE
            )
            if (affichage_output) {
              write.table(
                ecrit,
                file = stdout(),
                append = TRUE,
                quote = FALSE,
                sep = " ",
                na = "-",
                row.names = FALSE,
                col.names = FALSE
              )
            }
          }
        }
        
        else{
          # If the satellite is not visible from the receiver, due to the visibility masks
          ecrit <-
            t(c(
              ' The satellite',
              as.character(mon_data[num_ligne, 1]),
              'is not visible from the receiver.'
            ))
          write.table(
            ecrit,
            file = gsub(' ', '', paste(path_output, "output.txt")),
            append = TRUE,
            quote = FALSE,
            sep = " ",
            na = "-",
            row.names = FALSE,
            col.names = FALSE
          )
          if (affichage_output) {
            write.table(
              ecrit,
              file = stdout(),
              append = TRUE,
              quote = FALSE,
              sep = " ",
              na = "-",
              row.names = FALSE,
              col.names = FALSE
            )
          }
        }
        
        num_ligne <- num_ligne + 1
      }
    }
    
    close(progressbar) # We close the progress bar
    
    if (reflexion) {
      # If there were more than one reflection == if at least one satellite was visible and that the calculation converged
      
      # Creation of a text file containing the informations about the receiver --> 'receiver.txt'
      
      if (static_receiver) {
        pos_R_geo <- pos2geo(R_0_WGS84, R_p_WGS84, pos_R_xyz)
        write.table(
          'Longitude (o) - Latitude (o) - Ellipsoidal height (m) - vertical mask (o) - horizontal mask (o)',
          gsub(' ', '', paste(path_output, "receiver.txt")),
          append = FALSE,
          quote = FALSE,
          sep = " ",
          na = "-",
          row.names = FALSE,
          col.names = FALSE
        )
        write.table(
          t(c(
            round(pos_R_geo[1], 7),
            round(pos_R_geo[2], 7),
            round(pos_R_geo[3], 3),
            param[12, 2],
            param[13, 2],
            param[14, 2],
            param[15, 2]
          )),
          gsub(' ', '', paste(path_output, "receiver.txt")),
          append = TRUE,
          quote = FALSE,
          sep = " ",
          na = "-",
          row.names = FALSE,
          col.names = FALSE
        )
      }
      
      
      #
      # We transform all the text files in KML files in order to be importable with Google Earth
      #
      
      liste_output <- list.files(path_output)
      
      source('fenetre_info_kml.r') # We open a window to let the user know that the computer is computing...
      
      # We transform all the TXT files into KML files
      liste_output <- list.files(path_output)
      
      for (i in 1:length(liste_output)) {
        if (liste_output[[i]] == 'satellites_gps.txt' |
            liste_output[[i]] == 'satellites_glonass.txt' |
            liste_output[[i]] == 'satellites_galileo.txt' |
            liste_output[[i]] == 'receiver.txt' |
            liste_output[[i]] == 'specular_points_gps.txt' |
            liste_output[[i]] == 'specular_points_glonass.txt' |
            liste_output[[i]] == 'specular_points_galileo.txt' |
            liste_output[[i]] == 'traj_ondes_reflechies.tmp' |
            liste_output[[i]] == 'traj_ondes_directes.tmp')
          #if(liste_output[[i]]=='specular_points_gps.txt')
        {
          txt2kml('Simulation', liste_output[[i]], R_0_WGS84, R_p_WGS84)
        }
      }
      
      # We then concatenate all the different KML files (== we put the header, the styles, the main part, and the footer all together)
      source('concat_kml.r')
      liste_output <- list.files(path_output)
      concat_kml(path_output, "visible_satellite")
      concat_kml(path_output, "specular_points")
      concat_kml(path_output, "traj_ondes_reflechies")
      concat_kml(path_output, "traj_ondes_directes")
      concat_kml(path_output, "trace_sat")
      concat_kml(path_output, "all_sat")
      if (dessiner_fresnel == TRUE) {
        concat_kml(path_output, "fresnel_surface")
      }
      tkdestroy(ttinfo)
      
      tkdestroy(ttinfo_kml) # We close the information window
      
    }
    
    else {
      # If no satellite was visible, ever
      tkmessageBox(
        message = " No satellite was visible !",
        icon = "error",
        type = "ok",
        title = "Erreur !"
      )
      ecrit <- t(c('No satellite was visible !'))
      write.table(
        ecrit,
        file = gsub(' ', '', paste(path_output, "output.txt")),
        append = TRUE,
        quote = FALSE,
        sep = " ",
        na = "-",
        row.names = FALSE,
        col.names = FALSE
      )
      if (affichage_output) {
        write.table(
          ecrit,
          file = stdout(),
          append = TRUE,
          quote = FALSE,
          sep = " ",
          na = "-",
          row.names = FALSE,
          col.names = FALSE
        )
      }
    }
    
    system(
      'wine simu_del_temp.bat',
      show.output.on.console = FALSE,
      ignore.stderr = TRUE,
      ignore.stdout = TRUE
    ) # We delete all the temporary files we created
    
    # We display the total calculation time
    T2 <- Sys.time()
    print(difftime(T2, T1))
    #write.table(difftime(T2,T1), file = paste(path_output,'time_diff.txt',sep=""), append = TRUE, quote = FALSE, sep = " ", na = "-", row.names = FALSE, col.names = FALSE)
    
    # We warn the user that the calculations are finishe
    #tkmessageBox(message=" Simulation finished ! You can find all the output files in the folder Output/ ",icon="info",type="ok", title="Simulation finished !")
    
  }

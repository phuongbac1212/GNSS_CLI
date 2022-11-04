
### Calculation of E_sit (real elevation)

source('calc_amfx.r')
source('amfx_prlx.r')
nom_fichier_amf <- paste("AMF-X-F3A4C5G2-COR1_COR2_COR3_LEM1_LEM2-22921.3750_23130.5000.amf",sep='/')
station <- 'COR1'
date_jul_cnes <-  22921.5


E_sat <- 15.1614042784168*pi/180 # satellite elevation (rad)
A_sat <- 224.889657621929*pi/180# satellite azimuth (rad)
T_sat <- 94.6666582*pi/180 # satellite geocentric colatitude (rad)
L_sat <- -31.655568*pi/180 # satellite longitude  (rad)
R_sat <- 26270336.7562102  # satellite radius (geocentric coordinates)
T_sit <- 43.5929278*pi/180# geocentric colatitude of the site(rad)
F_sit <- 46.4070722*pi/180# geodetic latitude of the site (rad)
R_sit <- 6367392.48763499# site radius (geocentric coordinates)
L_sit <- 6.7167498*pi/180 # longitude of the site (rad)
R_sta <- 6367442.05052816
R_top <- 6450000

E_sit <- calc_prlx(nom_fichier_amf,station,date_jul_cnes,E_sat,A_sat,T_sat,L_sat,R_sat,T_sit,F_sit,R_sit,L_sit,R_top)[1]

nom_fichier_amf <- paste("AMF-X-F3A4C5G2-COR1_COR2_COR3_LEM1_LEM2-22921.3750_23130.5000.amf",sep='/')


A_sta <- A_sat + pi
R_top <- R_sta
D_top <- E_sit - calc_amfx(nom_fichier_amf,1,station,date_jul_cnes,E_sit,A_sta)*pi/180
P_top <- acos(R_sit/R_top*cos(D_top))-D_top
G_top <- sqrt(R_sit^2+R_top^2-2*R_sit*R_top*cos(P_top))



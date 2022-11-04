##################################################################
#
# janvier 2016 - Nicolas ROUSSEL
#
##################################################################
liste_codes_possibles <<- c('S1C','S1S','S1L','S1X','S1P','S1W','S1Y','S1M','S1N','S2C','S2D','S2S','S2L','S2X','S2P','S2W','S2Y','S2M','S2N','S5I','S5Q','S5X','S3I','S3Q','S3X','S1A','S1B','S1Z','S7I','S7Q','S7X','S8I','S8Q','S8X','S6A','S6B','S6C','S6Z') # liste des codes a rechercher dans le fichier d'observation
constante_gravitationnelle <<- 3.986004419*10^14      #constante gravitationnelle en m3/s2
vitesse_rotation_terre <<- 7.2921151467*10^-5    # Vitesse de rotation de la Terre en rad/sec
c_lumiere <<- 299792458
longueur_onde_L1_G <<- 0.1903 # en m
longueur_onde_L2_G <<- 0.2445  # en m
longueur_onde_L5_G <<- 0.2540  # en m
R_0_WGS84 <<- 6378137
R_p_WGS84 <<- 6356752.314245179497563967
f_WGS84 <<- 1/298.257223563 

#
# Determination du Work Directory et des differents chemins d'acces
#


wd <- getwd()
setwd('..')
path_main <<- paste(getwd(),'/',sep='')
path_input <<- paste(getwd(),'/0.input/',sep='')
path_extraction <<- paste(getwd(),'/1.extraction/',sep='')
path_output <<- paste(getwd(),'/2.output/',sep='')
path_Z_zip <<- paste(getwd(),'/programs/7-Zip/',sep='')
setwd(paste(getwd(),'/programs',sep='')) # We put back the Work Directory in the Programs folder


##################################################################
#
# Calculate the satellites coordinates from an almanac file, on a given date
#
# Launch: source("alm2xyz.R")
#
# 2nd semester 2012
#
# author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################



alm2xyz <- function(id,e,toa,i,OMEGA_point,racine_a,OMEGA_0,omega,M0,af0,af1,week_alm,GPS_semaine,jour,heure,minute,seconde,vitesse_rotation_terre,constante_gravitationnelle){

#
# id : PRN number of the satellite
# e : excentricity
# toe : Time reference of the ephemeris (seconds of the GPS week)
# i : orbit inclination
# OMEGA_point --> Rate of the right ascension of the node (rad / sec)
# racine_a : root of the half-major axis of the ellipsoid
# OMEGA_0 : Right ascension of the ascending node
# omega : argument of perigee
# M0 : mean anomaly at TOE
# af0 : Satellite clock bias
# af1 : Satellite clock drift
# week_alm : Number of the almanach week
# t0 : hour of the first measurement, taking the day number (sec) into account
# vitesse_rotation_terre : Speed of the Earth rotation (rad/sec)
# constante_gravitationnelle : GM (m3/s2)
#

theta <- as.numeric(vitesse_rotation_terre) # speed of the Earth rotation (rad/sec) --> cf http://hpiers.obspm.fr/eop-pc/index.php?index=constants&lang=en
GM <- as.numeric(constante_gravitationnelle) #m3s-2  (cf IERS Techn Notes > Tech Note 32 > Numerical Standard) Value of the gravitational constant of the earth
week_GPS <- week_alm + 1024 # 1024 is the gap between the alm and GPS weeks
t0 <- seconde + minute*60 + heure*60*60 + jour*24*60*60 + abs(GPS_semaine - week_GPS)*7*24*60*60
if (GPS_semaine < week_GPS) { t0 <- - t0}

a <- racine_a*racine_a
delta_t <- t0 - toa

# Determination of M

n0 <- sqrt(GM/a^3)
n <- n0 # + delta_n
M <- M0 + n*(delta_t)


# Determination of E
E <- M # We initialize E with M (radians !!!)
for (incr in 1:10){ # 10 iterations
    E_old <- E %% (2*pi)
    E <- (M + e*sin(E) ) %% (2*pi)
    dE <- (E-E_old) %% (2*pi)
    if (abs(dE)<10^-12){
        break
    }
}


# Determination of v
v <- atan2(sqrt(1-e^2)*sin(E),(cos(E)-e));


# ... some calculations ...
u0 <- omega + v
w <- omega #+ Cuc*cos(2*u0) + Cus*sin(2*u0)
u <- w + v
r0 <- a*(1-e*cos(E))
r <- r0 #+ Crc*cos(2*u0) + Crs*sin(2*u0)

# Calculating Cartesian coordinates related to the orbital plane
x_prime <- r*cos(u)
y_prime <- r*sin(u)
z_prime <- 0


# Determination of i and OMEGA
#i <- io + Cic*cos(2*u0) + Cis*sin(2*u0) + idot*(delta_t)
OMEGA <- OMEGA_0 + (OMEGA_point - theta)*(delta_t) - theta * toa

# Calculation of Cartesian geocentric coordinates

X <- x_prime*cos(OMEGA) - y_prime*cos(i)*sin(OMEGA)
Y <- x_prime*sin(OMEGA)+ y_prime*cos(i)*cos(OMEGA)
Z <- y_prime*sin(i)

return(c(X,Y,Z)/1000) # in km
}

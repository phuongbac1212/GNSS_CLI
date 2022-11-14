##################################################################
#
# Function to determine the size of the footprint of the GNSS signal (1st Fresnel surface)
#
# Launch: source("surf_Fresnel.R")
#
# 2nd semester 2012
#
# Author : Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################

surf_Fresnel <- function(lambda, hauteur, elev){
    # lambda: wave length (m)
    # hauteur: height of the receiver (m)
    # elev : satellite elevation (rad)
    lambda <- as.numeric(lambda)
    hauteur <- as.numeric(hauteur)
    elev <- as.numeric(elev)
    
    axes_Fresnel <- c(0,0,0)
    
    #
    # Calculation of the semi-minor and semi-major axis and the area of the first Fresnel surface
    #
      if(hauteur>=2000){
    axes_Fresnel[1] <-  sqrt(lambda * hauteur * sin(elev)) / sin(elev) ^ 2
    axes_Fresnel[2] <-  sqrt(lambda * hauteur * sin(elev)) / sin(elev)
    axes_Fresnel[3] <- pi * axes_Fresnel[1] * axes_Fresnel[2]
        write.table(axes_Fresnel[1],'a.txt',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(axes_Fresnel[2],'b.txt',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
       }
    else{
    d <- lambda/2
    b <- sqrt((2*d*hauteur)/sin(elev)+(d/sin(elev))^2)
    a <- b/sin(elev)
    write.table(a,'a.txt',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(b,'b.txt',append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
    axes_Fresnel[1] <- a
    axes_Fresnel[2] <- b
    axes_Fresnel[3] <-  pi * axes_Fresnel[1] * axes_Fresnel[2]
     }
    return(axes_Fresnel)

}
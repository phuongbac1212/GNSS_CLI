
fresnel_filter <- function(x) {
  nb_gones <- 40 # How many sides of the ellipse (= approximated as a polygone)
  nb_gones <- nb_gones +2
  # Coordinates of the known points in the 3D WGS84 system
  
  a = x$Fz_a
  b = x$Fz_b
  
  source("src/func/geo2pos.R")
  pos_S_xyz = geo2pos(cbind(x$Lon_simu, x$Lat_simu, x$Height_simu))
  
  temp = c(x$A_x,x$A_y, x$A_z) # cal cooord_A
  x1 <- temp[1] # coord_A
  y1 <- temp[2]
  z1 <- temp[3]
  
  temp = c(x$B_x,x$B_y, x$B_z) # cal cooord_b
  x2 <- temp[1] # coord_B
  y2 <- temp[2]
  z2 <- temp[3]
  
  temp = c(x$AB_x,x$AB_y, x$AB_z) # cal cooord_aa_bb
  x3 <- temp[1] # coord_AA_BB
  y3 <- temp[2]
  z3 <- temp[3]

  x0 <- pos_S_xyz[1]  # coord_S
  y0 <- pos_S_xyz[2]
  z0 <- pos_S_xyz[2]
  
  # Coordinates of the known points in the 3D local system, defined by X : axis S->A, Y : axis S->B, Z : the straight line orthogonal to forme a direct trihedral.
  
  X1 <- a
  Y1 <- 0
  Z1 <- 0
  
  X2 <- 0
  Y2 <- b
  Z2 <- 0
  
  X3 <- -a
  Y3 <- 0
  Z3 <- -b
  
  X0 <- 0
  Y0 <- 0
  Z0 <- 0
  
  u <- 0
  v <- 0
  
  
  
  # We look for the parameters to go from a system to the other.
  # To do that, we do a 9 parameters transformation :
  # We have : 
  # x = T1 + A*X + B*Y + C*Z
  # y = T2 + D*X + E*Y + F*Z
  # z = T3 + G*X + H*Y + I*Z
  # and so we have: A %*% X = B, with :
  # A =   X1 Y1 Z1 0 0 0 0 0 0
  #       0 0 0 X1 Y1 Z1 0 0 0
  #       0 0 0 0 0 0 X1 Y1 Z1
  #       X2 Y2 Z2 0 0 0 0 0 0
  #       0 0 0 X2 Y2 Z2 0 0 0
  #       0 0 0 0 0 0 X2 Y2 Z2
  #       X3 Y3 Z3 0 0 0 0 0 0
  #       0 0 0 X3 Y3 Z3 0 0 0
  #       0 0 0 0 0 0 X3 Y3 Z3
  # X =  c(A,B,C,D,E,F,G,H,I)
  # B = (x1,y1,z1,x2,y2,z2,x3,y3,z3)
  
  t1 <- pos_S_xyz[1]
  t2 <- pos_S_xyz[2]
  t3 <- pos_S_xyz[3]
  
  A <- matrix(c(X1,Y1,Z1,0,0,0,0,0,0,0,0,0,X1,Y1,Z1,0,0,0,0,0,0,0,0,0,X1,Y1,Z1,X2,Y2,Z2,0,0,0,0,0,0,0,0,0,X2,Y2,Z2,0,0,0,0,0,0,0,0,0,X2,Y2,Z2,X3,Y3,Z3,0,0,0,0,0,0,0,0,0,X3,Y3,Z3,0,0,0,0,0,0,0,0,0,X3,Y3,Z3),ncol=9,byrow=TRUE)
  B <- c(x1-t1,y1-t2,z1-t3,x2-t1,y2-t2,z2-t3,x3-t1,y3-t2,z3-t3) 
  
  result <- qr.solve(A,B)
  A <- result[1]
  B <- result[2]
  C <- result[3]
  D <- result[4]
  E <- result[5]
  F <- result[6]
  G <- result[7]
  H  <- result[8]
  I <- result[9]   
  
  # XX : abscissa of the differents points of the ellipse we will plot
  XX1 <- seq(from=u-a, to =u+a, length =round(nb_gones/2))
  XX2 <- seq(from=u+a, to =u-a, length =round(nb_gones/2))
  XX <- c(XX1,XX2)
  # YY : ordinates of the differents points of the ellipse we will plot
  YY1 <- sqrt(b^2-((XX1-u)^2/a^2*b^2)) + v
  YY2 <- -sqrt(b^2-((XX2-u)^2/a^2*b^2)) + v
  YY <- c(YY1,YY2)
  ZZ <- seq(from=0,to=0,length=length(YY))
  xx <- 0
  yy <- 0
  zz <- 0
  long <- 0
  lat <- 0
  he <- 0
  # We transform the points coordinates of the ellipse we will plot from the 3D local system to the 3D WGS84 system thanks to the 9 parameters we determined before.
  for(u in 1:length(XX)){
    xx[u] <- t1 + A*XX[u] + B*YY[u] + C*ZZ[u]
    yy[u] <- t2 + D*XX[u] + E*YY[u] + F*ZZ[u]
    zz[u] <- t3 + G*XX[u] + H*YY[u] + I*ZZ[u]   
    result <- pos2geo(c(xx[u],yy[u],zz[u]))   
    long[u] <- result[1]
    lat[u] <- result[2]
    he[u]<- result[3]
  }
  if(x2==y2 & y2==z2 & z2==0){ # If an error occured during the calculation of B (sometimes, the discriminant is negative, I don't know why) --> we don't plot this Fresnel surface
    long<- x0
    lat <- y0
  }
  
  shape.file = "roi"
  pol = read_sf(dsn = shp_path, layer = shape.file)
  pol.coor = as.data.frame(st_coordinates(pol)[,c("X", "Y")])
  
  point = point.in.polygon(long, lat, pol.coor$X, pol.coor$Y)
  return (sum(which(point == 0)) == 0)
}
# # test
# tmp = NA
# for (i in 1:nrow(res)) {
#   print(res[i,])
#   tmp[i] <- fresnel_filter(res[i,])
# }
# res = res[tmp,]

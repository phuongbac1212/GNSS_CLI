require(splus2R)

calcul_normale <- function(pt1,pt2){
  
  xA <- pt1[1]
  yA <- pt1[2]
  zA <- pt1[3]
  xB <- pt2[1]
  yB <- pt2[2]
  zB <- pt2[3]
  
  if(xA==xB){xB<-xB+0.0000000001}
  if(yA==yB){yB<-yB+0.0000000001}
  if(zA==zB){zB<-zB+0.0000000001}
  # plane going through the 2 points and the center of the Earth
  A0 <- c(-xA,-yA,-zA)
  
  AB <- c(xB-xA,yB-yA,zB-zA)
  
  vecto <- c(A0[2]*AB[3]-A0[3]*AB[2],A0[3]*AB[1]-A0[1]*AB[3],A0[1]*AB[2]-A0[2]*AB[1])
  
  a <- vecto[1]
  b <- vecto[2]
  cc <- vecto[3]
  d <- -xA*vecto[1]-yA*vecto[2]-zA*vecto[3]
  # plane normal to the straight line (pt1 pt2)
  A <- (xB-xA)
  B <- yB-yA
  CC <- zB-zA
  DD <- -xA*A-yA*B-zA*CC
  
  # intersection of the two planes
  x <- (-d+DD*b/B)/(a-A*b/B)
  y <- (-DD-A*x)/B
  
  normale <- c((x-xA),(y-yA),0-zA)
  normale <- normale/vecnorm(normale)
  return(normale)
}

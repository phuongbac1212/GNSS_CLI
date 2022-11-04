##################################################################
#
# Interpolateur de Lagrange
#
# lancement: source("post_lagrange.R")
#
# 2eme Semestre 2012
#
##################################################################
 
lagrange <- function(x,y,x0){

#
# x : vecteur contenant les abscisses des points connus
# y : vecteur contenant les ordonnees des points connus
# x0 : vecteur contenant les abscisses des points dont on cherche l'ordonnee
#
  
library(polynom)
x <- as.numeric(x)
y <- as.numeric(y)
print(poly.calc(as.numeric(x[1:100]),as.numeric(y[1:100]),tol=100))
zef
return(predict(poly.calc(as.numeric(x),as.numeric(y),tol=10^-2),x0))

}
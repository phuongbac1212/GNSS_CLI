intersection_droite_ellipsoide <- function(a,b,A,B,C,xa,ya,za){
K1 <- A^2/a^2+B^2/a^2+C^2/b^2
K2 <- 2*A*xa/a^2+2*B*ya/a^2+2*C*za/b^2
K3 <- xa^2/a^2+ya^2/a^2+za^2/b^2-1

delta <- K2^2-4*K1*K3
t1 <- (-K2+sqrt(delta))/(2*K1)
t2 <- (-K2-sqrt(delta))/(2*K1)

x1 <- A*t1+xa
y1 <- B*t1+ya
z1 <- C*t1+za
x2 <- A*t2+xa
y2 <- B*t2+ya
z2 <- C*t2+za

d1 <- vecnorm(c(x1,y1,z1)-c(xa,ya,za))
d2 <- vecnorm(c(x2,y2,z2)-c(xa,ya,za))

if(d1<d2){
res <- c(x1,y1,z1)
}
else{
res <- c(x2,y2,z2)
}
return(res)
}
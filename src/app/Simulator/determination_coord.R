##################################################################
#
# Function to determine the point B (Xb,Yb,Zb) :
#   1. belonging to a plane (ASB) :  uXb + vYb + wZb + r = 0   *
#   2. belonging to a sphere centered in S (Xo,Yo,Zo) and with a radius b      (Xb-Xo)� + (Yb-Yo)� + (Zb-Zo)� = b�  
#   3. SA.SB = 0 ==> (Xb - Xo)(Xa - Xo) + (Yb-Yo)(Ya -Yo) + (Zb-Zo)(Za - Zo) = 0
#
# Launch: source("determination_coord.R")
#
# 2nd semester
#
# Author: Nicolas ROUSSEL - nicolas.roussel@get.obs-mip.fr
#
##################################################################


determination_coord <- function(Xa,Ya,Za,u,v,w,r,Xo,Yo,Zo,b){
  # print(paste("W: ", Xa,Ya,Za,u,v,w,r,Xo,Yo,Zo,b))
  # print(paste("type", typeof(Xa),typeof(Ya),typeof(Za),typeof(u),typeof(v),typeof(w),typeof(r),typeof(Xo),typeof(Yo),typeof(Zo),typeof(b)))
  options(digits=22)
  
  # We wants to find the point belonging to :
  #    - the local plane uXb + vYb + wZb + r = 0 (2)
  #    - the sphere centered in S and with a radius b : (Xb-Xo)� + (Yb-Yo)� + (Zb-Zo)� = b�  (3)
  # Moreover, the scalar product SA.SB = 0, so : (Xb - Xo)(Xa - Xo) + (Yb-Yo)(Ya -Yo) + (Zb-Zo)(Za - Zo) = 0 (1)
  
  
  
  #           (2)
  #           Xb = -vYb/u - wZb/u - r/u
  #          
  #          (2) --> (1)        
  #           (-vYb/u - wZb/u - r/u - Xo)(Xa - Xo) + (Yb-Yo)(Ya -Yo) + (Zb-Zo)(Za - Zo) = 0    (11)
  #          
  #           (2) --> (3)          
  #           (-vYb/u - wZb/u - r/u - Xo)� + (Yb-Yo)� + (Zb-Zo)� = b�      (33)
  #            
  #          
  #          (11)         
  #            -vYbXa/u - wZbXa/u - rXa/u - XoXa + vYbXo/u + wZbXo/u +rXo/u + Xo� + YbYa-YbYo - YoYa + Yo� + (Zb-Zo)(Za-Zo) = 0 
  #            ==>           
  #            Yb = [wZbXa/u + rXa/u + XoXa - wZbXo/u - rXo/u -Xo� +YoYa - Yo� -(Zb-Zo)(Za-Zo)]/[-vXa/u+vXo/u+Ya - Yo]
  #          
  #          (11) --> (33)                                                           
  #            (([-vwZbXa/u� - vrXa/u� - vXoXa/u + vwZbXo/u� + vrXo/u� +vXo�/u -vYoYa/u + vYo�/u +vZbZa/u-vZbZo/u-vZoZa/u+vZo�/u]/[-vXa/u+vXo/u+Ya - Yo])  - wZb/u - r/u - Xo)�
  #            + (([wZbXa/u + rXa/u + XoXa - wZbXo/u - rXo/u -Xo� +YoYa - Yo� - ZbZa+ZbZo+ZoZa-Zo�]/[-vXa/u+vXo/u+Ya - Yo])-Yo)� + Zb� -2ZbZo + Zo�= b�                                                                                                                                                       
  #          
  #              (A1*Zb + K1)� + (A2*Zb + K2)�   +A3*Zb + K3 + Zb�= 0
  
  #             A1 = (-vwXa/u�+vwXo/u�+vZa/u-vZo/u)/(-vXa/u + vXo/u + Ya - Yo) -w/u
  #            K1 = (-vrXa/u�-vXoXa/u+vrXo/u�+vXo�/u-vYoYa/u+vYo�/u-vZoZa/u+vZo�/u)/(-vXa/u + vXo/u + Ya - Yo) - r/u - Xo
  #            A2 = (wXa/u -wXo/u -Za +Zo)/(-vXa/u+vXo/u+Ya - Yo)
  #            K2 = (rXa/u + XoXa - rXo/u - Xo� + YoYa - Yo� +ZoZa - Zo�)/(-vXa/u+vXo/u+Ya - Yo)  -Yo
  #            A3 = -2Zo
  #            K3 = Zo� - b�
  #            
  #          (A1*Zb + K1)� + (A2*Zb + K2)�   +A3*Zb + K3 + Zb�= 0
  #              ==> (A1�Zb� + 2K1A1Zb + K1� + A2�Zb� + K2� + 2A2ZbK2 + A3Zb + K3 + Zb� = 0
  #           
  #           ==> B1*Zb�  +  B2*Zb  +  B3 = 0
  #           
  #           B1 =  A1� + A2� + 1
  #           B2 = 2K1A1 + 2A2K2 + A3
  #            B3 = K1� + K2� + K3 
  
  A1 <- (-v*w*Xa/(u*u)+v*w*Xo/(u*u)+v*Za/u-v*Zo/u)/(-v*Xa/u + v*Xo/u + Ya - Yo) -w/u
  K1 <- (-v*r*Xa/(u*u)-v*Xo*Xa/u+v*r*Xo/(u*u)+v*(Xo*Xo)/u-v*Yo*Ya/u+v*(Yo*Yo)/u-v*Zo*Za/u+v*Zo*Zo/u)/(-v*Xa/u + v*Xo/u + Ya - Yo) - r/u - Xo
  A2 <- (w*Xa/u -w*Xo/u -Za +Zo)/(-v*Xa/u+v*Xo/u+Ya - Yo)
  K2 <- (r*Xa/u + Xo*Xa - r*Xo/u - Xo*Xo + Yo*Ya - Yo*Yo +Zo*Za - Zo*Zo)/(-v*Xa/u+v*Xo/u+Ya - Yo)  -Yo
  A3 <- -2*Zo
  K3 <- Zo*Zo - b*b
  
  B1 <- A1*A1 + A2*A2 + 1
  B2 <- 2*K1*A1 + 2*A2*K2 + A3
  B3 <- K1*K1 + K2*K2 + K3
  
  # Resolution of the equation B1*Zb�  +  B2*Zb  +  B3 = 0
  
  discriminant <- B2*B2 - 4*B1*B3
  
  if(discriminant <0){
    print('Error during the claculation of the first Fresnel surface')
    B <- c(0,0,0)
    BB <- c(0,0,0) 
    discriminant <- -discriminant
  }
  
  else{
    Zb1 <- (-B2 - sqrt(discriminant)) / (2*B1)
    Zb2 <- (-B2 + sqrt(discriminant)) / (2*B1)
    Yb1 <- (w*Zb1*Xa/u + r*Xa/u + Xo*Xa - w*Zb1*Xo/u - r*Xo/u -Xo^2 +Yo*Ya - Yo^2 -(Zb1-Zo)*(Za-Zo))/(-v*Xa/u+v*Xo/u+Ya - Yo)
    Yb2 <- (w*Zb2*Xa/u + r*Xa/u + Xo*Xa - w*Zb2*Xo/u - r*Xo/u -Xo^2 +Yo*Ya - Yo^2 -(Zb2-Zo)*(Za-Zo))/(-v*Xa/u+v*Xo/u+Ya - Yo)
    Xb1 <- -v*Yb1/u - w*Zb1/u - r/u
    Xb2 <- -v*Yb2/u - w*Zb2/u - r/u
    B <- c(Xb1,Yb1,Zb1)
    BB <- c(Xb2,Yb2,Zb2)
  }
  
  return(c(B,BB))
}
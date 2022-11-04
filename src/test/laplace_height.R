#localPb = local pressure in mB
#heightm = height in m
#tempC = temparature in Cencilous deg 
#wpppmB =  Station water vapour pressure in mB
#latit = latitude in deg
laplace.correction = function(localPb ,
                          heightm,
                          tempC = 0,
                          wppmB = 0,
                          latit = 0) {
  # Laplace
  p = localPb
  z = heightm
  p0 = p * 10 ^ (z / 18400)
  
  t = tempC
  pt = p * 10 ^ (z / (18400 * (1 + (t * 0.00370))))
  hz = wppmB
  hf = 1 / (1 - 0.378 * (hz / p))
  phz = p * 10 ^ (z / (18400 * (1 + (t * 0.00370)) * hf))
  
  lz  = latit
  hl = 1 / (1 - 0.0026 * (cos(2 * lz / 57)))
  plz = p * 10 ^ (z / (18400 * (1 + (t * 0.00370)) * hf * hl))
  
  gz	= 1 + z / (2 * 6367324)
  pgz = p * 10 ^ (z / (18400 * (1 + (t * 0.00370)) * hf * hl * gz))
  
  # ##BABINET'S FORMULA
  # a = (16000 + (64 * t))
  # az = a * 1 + z * 1
  # am = a * 1 - z * 1
  # ax = az / am
  # 
  # Pbab = (p * (ax))
  
  return(pgz)
}


hypsometric = function(P0, P, t) {
  return(((P0/P)^(1/5.257)-1)*(t+273.15)/0.0065)
}

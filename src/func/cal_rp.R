cal.frzone <- function(cx, cy, a, b, rotAngle) {
  # # This functions returns points to draw an ellipse
  # %#
  # %#  @param cx     center x position of the reflected point
  # %#  @param cy     center y position of the reflected point
  # %#  @param a     Semimajor axis
  # %#  @param b     Semiminor axis
  # %#  @param angle Azimuth of the FZ (in degrees)

  step = 30
  angle = pracma::linspace(0, 2*pi, step)
  
  X = a * cos(angle);
  Y = b * sin(angle);
  xRot = X*cos(rotAngle) - Y*sin(rotAngle);
  yRot = X*sin(rotAngle) + Y*cos(rotAngle);
  
  return(list(cbind("X" = cx+xRot, "Y" = cy+yRot)))
}

cal.rp.pos <- function(lat, lon, h, h0, f, lambda, elev, azi) {
  h.m = h0 - (lambda*f)/2 
  d = lambda/2
  b = sqrt((2*d*h)/sin(elev)+(d/sin(elev))^2) 
  a = b/sin(elev)
  
  dist.rp = h0/azi
  dx = dist.rp*sin(azi)
  dy = dist.rp*cos(azi)
  dlon = dx/(111320*cos(lat))
  dlat = dy/110540
  
  lat.rp = lat + dlat
  lon.rp = lon + dlon
  h.rp = h - h0
  
  source("src/func/geo2pos.R")
  source("src/func/pos2geo.R")
  rp.xyz = geo2pos(cbind(lon.rp, lat.rp, h.rp))
  fz = cal.frzone(rp.xyz[1], rp.xyz[2], a, b, azi)
  fz[[1]] = cbind(fz[[1]], "Z" = rp.xyz[3])
  return(fz)
}


#for test
require(fst)
require(dplyr)
require(data.table)
x = data.frame(fst("output/extract/backup/2205_2_PC06_0.fst"))
rec.ll = pos2geo(rec)
x$f = 299792458/ as.numeric(x$Wavelength)

test  = apply(
  x[1:100,],
  1,
  FUN = function(t) {
    print(t)
    a = cal.rp.pos(
      rec.ll[2],
      rec.ll[1],
      rec.ll[3],
      h0,
      as.numeric(t["f"]),
      as.numeric(t["Wavelength"]),
      as.numeric(t["elevation"]),
      as.numeric(t["azimuth"])
    )
    b = t(as.data.frame(apply(a[[1]], 1, pos2geo)))
    return(list(b))
  }
)

# require(leaflet)
# leaflet() %>% addTiles() %>% addPolygons(lng = b[,1], lat=b[,2])
require(ggplot2)
plot = ggplot() +
  geom_point( aes(x=rec.ll[1], y = rec.ll[2]))

for (i in 1:length(test)) {
  geo = test[[i]][[1]]
  plot = plot + geom_polygon(aes(geo[,1], geo[,2]))
}
plot

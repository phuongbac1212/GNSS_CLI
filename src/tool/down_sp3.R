source("src/lib/sp3_handler.R")
sp3.url = "ftp://igs.ign.fr/pub/igs/products/mgex/"
down_sp3 <- function(from, to) {
  
  d.from  = date(as_datetime(from, origin = "1980-01-06"))
  d.to =  date(as_datetime(to, origin = "1980-01-06"))
  while (d.from < d.to) {
    link = sp3.getlink(as.numeric(as.POSIXct(d.from, origin = "+0")) - 315964782)
    sp3.download(link, "SP3/")
    d.from = d.from +1
  }
}
#wtc
my.outlier <- function(x, c.name, coef, bias) {
  # x = x[, c.name]
  x = as.data.frame(x)
  c.name = match(c.name, names(x))
  if (!is.na(coef)) {
    out = x[, c.name[2]] %in% boxplot.stats(x[, c.name[2]], coef = coef)$out
    x = x[!out,]
  }
  if (!is.na(bias)) {
    spline =  lm(x[, c.name[2]] ~ splines::bs(x[, c.name[1]], 100))$fitted.values
    
    diff = abs(spline - x[, c.name[2]])
    out2 = diff <= bias
    x = x[out2,]
  }
  return(x)
}

require(Metrics)
require(tsutils)
require(ggplot2)
require(lubridate)
require(tidyr)
require(dplyr)
require(Rwave)
require(biwavelet)
res = read.csv("output/result/Result.csv") %>% drop_na(height)
res$time = lubridate::as_datetime(res$timestamp, origin = "1980-1-6")
res$height = 7.5 - res$height
res = filter(res, height != max(height))
res = my.outlier(res, c("time", "height"), 2, 2) %>% drop_na(height)
plot(res$time, res$height, type = "p")
res$group = round(res$time, "hour")

gnss.data = res[, c("height", "time", "group")] %>% group_by(group) %>% 
  summarise_at(vars(height), list(h_mean_hour = mean))
gnss.data$time = as.POSIXct(gnss.data$group)

list.me007 = list.files("input/sensor/", full.names = T)
sensor2 = lapply(list.me007, read.csv, header=F) %>% data.table::rbindlist()
names(sensor2) <- c("time", "height")
sensor2$time = lubridate::as_datetime(sensor2$time)

sensor2 = dplyr::filter(sensor2, height > 4000)
sensor2= my.outlier(sensor2, c("time", "height"), 2, 200)
sensor2 =  sensor2 %>% group_by(round(time, "hour")) %>% 
  summarise_at(vars(height), list(height_mean_hour = mean))
names(sensor2) <- c("time", "height")
sensor2$height = 4.5 - sensor2$height/ 1000


a = approx(gnss.data, n= sum(diff(gnss.data$time))+1)
t = as.POSIXct(a$x, origin= "1970-1-1")
d1 = cbind(t, a$y)

a2 = approx(sensor2, n= sum(diff(sensor2$time))+1)
t2 =  as.POSIXct(a2$x, origin= "1970-1-1")
d2 = cbind(t2,a2$y)

d1 = as.data.frame(d1)
d2= as.data.frame(d2)
names(d1) = c("time", "height")
names(d2) = c("time", "height")
d3 = inner_join(d1, d2, by = "time")

xa = as.POSIXct(d1$V1, origin="1970-1-1")
x = xwt(d1, d2)
plot.biwavelet(x,plot.phase	= T, xaxt="n", yaxt="n")
xlocs <- pretty_dates(xa, 10)
axis(side = 1, at = xlocs, labels = format(xlocs, "%d/%m/%y"))
axis.locs <- axTicks(2)
yticklab <- format(round(2^axis.locs/60/60/24, 2))
axis(2, at = axis.locs, labels = yticklab)



cwt(d1, noctave = 10)
plot(wt(d1))

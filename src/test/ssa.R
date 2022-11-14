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

normalize <- function(arr, x, y) {
  m = min(arr)
  range = max(arr) - m
  arr = (arr - m)/range
  
  range2 = y-x
  norm= (arr *range2) +x
}
require(tsutils)
require(ggplot2)
require(lubridate)
require(tidyr)
require(dplyr)

res = read.csv("output/result/Result.csv") %>% drop_na(height)
res$time = lubridate::as_datetime(res$timestamp, origin = "1980-1-6")
res$height = 7.5 - res$height
res = filter(res, height != max(height))


########## location filter #############
library(leaflet)
library(randomcoloR)
p = unique(res$PRN)
col = randomColor(count = length(p))
a = data.frame(1:length(res$PRN), res$PRN)
names(a) = c("A", "B")
b = data.frame(col, p)
names(b) = c("A", "B")
a = inner_join(a, b, by = "B")
leaflet() %>% addProviderTiles(providers$Esri.WorldImagery) %>% addCircleMarkers(
  lng = res$Lon_simu,
  lat = res$Lat_simu,
  color = a$A.y,
  radius = 1
)

source("src/func/location_filter.R")
shp_path = "input/shp/"
res = location.filter(point.dataframe = res)



res = my.outlier(res, c("time", "height"), 2, 2) %>% drop_na(height)
res$group = substr(res$PRN, 2,2)

gnss.data = res[, c("height", "time", "group")] %>% group_by(group = round(res$time, "hour")) %>%
  summarise_at(vars(height), list(h_mean_hour = mean))
gnss.data$time = as.POSIXct(gnss.data$group)


# SSA
X = gnss.data$h_mean_hour
M = 14
N = length(X)
lag = lagmatrix(X, 0:-(M - 1))
lag = apply(lag, 2, replace_na, 0)
C = t(lag) %*% lag / 20
e = eigen(C)
lambda = e$values
RHO = e$vectors
PC = lag %*% RHO

library(pracma)
RC = zeros(N, M)
for (m in (1:M)) {
  Z = lagmatrix(PC[,m], 0:-(M - 1))
  Z = apply(Z, 2, replace_na, 0)
  RC[,m] = Z %*% RHO[,m] / M
}
################################################################################list.me007 = list.files("input/sensor/", full.names = T)
sensor2 = lapply(list.me007, read.csv, header = F) %>% data.table::rbindlist()
names(sensor2) <- c("time", "height")
sensor2$time = lubridate::as_datetime(sensor2$time, origin="1970-01-01")

sensor2 = dplyr::filter(sensor2, height > 4000)
sensor2 = my.outlier(sensor2, c("time", "height"), 2, 200)
sensor2 =  sensor2 %>% group_by(round(time, "hour")) %>%
  summarise_at(vars(height), list(height_mean_hour = mean))
names(sensor2) <- c("time", "height")
sensor2$height = 4.3 - sensor2$height / 1000

list.bme = list.files("input/bme/", full.names = T)
sensor = lapply(list.bme, read.csv, header = F) %>% data.table::rbindlist()
names(sensor) <- c("ts", "temprature", "pressure", "humid")

sensor$time = lubridate::as_datetime(sensor$ts, origin="1970-01-01")

sensor3 =  sensor %>% group_by(round(time, "hour")) %>%
  summarise_at(vars(humid), list(atm_mean_hour = mean))
names(sensor3) = c("time", "humid")

sensor4 =  sensor %>% group_by(round(time, "hour")) %>%
  summarise_at(vars(temprature), list(atm_mean_hour = mean))
names(sensor4) = c("time", "temprature")

sensor =  sensor %>% group_by(round(time, "hour")) %>%
  summarise_at(vars(pressure), list(atm_mean_hour = mean))
names(sensor) = c("time", "atm")

h = -(sensor$atm - mean(sensor$atm)) / (1.020 * 9.81)
h = data.frame(sensor$time, h)
################################################################################
require(easyGgplot2)
sensor = dplyr::filter(sensor, time >= min(gnss.data$time) &
                         time <= max(gnss.data$time))
g = ggplot() +
  theme(axis.text = element_text(size = 12), axis.title.x=element_blank()) +
  scale_x_datetime(breaks = seq(
    from = min(as.POSIXct(t)),
    to = max(as.POSIXct(t)),
    by = 3 * 24 * 60 * 60
  ),
  date_labels = "%d-%m") +
  scale_y_continuous(
    labels = function(x)
      sprintf("%10.1f", x)
  )


ggplot2.multiplot(
  g + geom_line(aes(y = X, x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 1], x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 2], x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 3], x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 4], x = gnss.data$time)),
  g + geom_line(aes(
    y = sensor2$height, x = as.POSIXct(sensor2$time)
  )),
  g + geom_line(aes(
    y = sensor$atm, x = as.POSIXct(sensor$time)
  )),
  cols = 1
)
gnss.data.old = gnss.data
ggplot2.multiplot(
  g + geom_line(aes(y = X, x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 1], x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 2], x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 3], x = gnss.data$time)),
  g + geom_line(aes(y = RC[, 4], x = gnss.data$time)),
  g + geom_line(aes(
    y = sensor2$height, x = as.POSIXct(sensor2$time)
  )) + geom_line(aes(x= gnss.data$time, y= gnss.no.pc1, col = "no pc1")),
  g + geom_line(aes(
    y = sensor$atm, x = as.POSIXct(sensor$time)
  )),
  cols = 1
)
##########################plot them             ################################
################################################################################
library(imputeTS)
library(xts)
library(plotly)
t_lim = c(min(gnss.data$time), max(gnss.data$time))


gnss.xts = xts(x = X, order.by = gnss.data$time)
sensor2 = filter(sensor2, as.numeric(time) > 1663365600 | as.numeric(time) <1662343200)

# m = inner_join(h, gnss.data, by = c("sensor.time" = "time"))
# cor(m$h_mean_hour, m$h)
# rmse(m$h, m$h_mean_hour)
# me007.xts = as.xts(sensor2$height, order.by = sensor2$time)
# plot.xts(me007.xts, type = "p")
# plot(me007.xts, type = "b")
# dygraph(me007.xts)

par(mfrow = c(6, 1),mar = c(2, 4, 0.5, 2), family= "serif", cex = 0.9) # puts 4 plots in one window (2x2)

library(Metrics)
plot(
  y = X,
  x = gnss.data$time,
  type = "l",
  xlab = "date",
  ylab = "Water height(m)",
  xlim = t_lim,
  xaxt = 'n',
  col = "red"
)
lines(sensor2$height, x = as.POSIXct(sensor2$time),  
     col ="darkblue")
xlocs <- pretty_dates(xa, 15)
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))
legend("topleft", legend=c("GNSS-based", "Sensor"),
       col=c("red", "darkblue"), lty=1:2, cex=0.8, bty = "n")


plot(y = PC[, 1], x = gnss.data$time,  
    type = "l",
  ylab = "PC1", xlim = t_lim, xaxt = 'n', col = "red")
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))

plot(y = PC[, 2], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC2", xlim = t_lim, xaxt = 'n', col = "red")
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))
plot(y = PC[, 3], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC3", xlim = t_lim, xaxt = 'n', col = "red")
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))
plot(y = PC[, 4], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))

plot(sensor$atm, x = as.POSIXct(sensor$time),  
    xlab = "Height",
     ylab = "Atmospheric(mbar)", xlim = t_lim, xaxt = 'n', type = "l", col = "red")
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))

################################################################################
################################################################################
################################################################################
plot(y = PC[, 5], x = gnss.data$time,  
type = "l", xlab = "Height",
ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 6], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 7], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 8], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 9], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 10], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 11], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 12], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 13], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
plot(y = PC[, 14], x = gnss.data$time,  
     type = "l", xlab = "Height",
     ylab = "PC4", xlim = t_lim, xaxt = 'n', col = "red")
################################################################################

# xwt
gnss.data = gnss.data.old
library(biwavelet)
library(lubridate)
#gnss.data$h_mean_hour =
a = approx(gnss.data, n = sum(diff(gnss.data$time)) + 1)
t = as.POSIXct(a$x, origin = "1970-1-1")
d1 = cbind(t, a$y)

a1 = approx(sensor2, n= sum(diff(sensor2$time))+1)
t1 = as.POSIXct(a1$x, origin = "1970-1-1")
d1 = cbind(t1, a1$y)

a2 = approx(h, n = sum(diff(h$sensor.time)) + 1)
t2 =  as.POSIXct(a2$x, origin = "1970-1-1")
d2 = cbind(t2, a2$y)

d1 = as.data.frame(d1)
d2 = as.data.frame(d2)
names(d1) = c("time", "height")
names(d2) = c("time", "height")
d3 = inner_join(d1, d2, by = "time")
d1 = data.frame(d3$time, d3$height.x)
d2 = data.frame(d3$time, d3$height.y)
names(d1) = c("time", "height")
names(d2) = c("time", "height")
xa = as.POSIXct(d1$time, origin = "1970-1-1")
x = xwt(d1, d2)
plot.biwavelet(x,
               plot.phase	= T,
               xaxt = "n",
               yaxt = "n")
xlocs <- pretty_dates(xa, 15)
axis(side = 1,
     at = xlocs,
     labels = format(xlocs, "%d/%m"))
yticklab = c("0.125", "0.25", "0.5", "1", "2", '4', '8', '16', '32')
axis.locs = log2(as.numeric(yticklab) * 24 * 60 * 60)
axis(2, at = axis.locs, labels = yticklab)


ggplot() +
  geom_line(size = 2, aes(y = a$y, x = as.POSIXct(t), colour = "GNSS-based")) +
  geom_line(size = 2, aes(
    y = a2$y,
    x = as.POSIXct(t2),
    colour = "Atmospheric"
  )) +
  xlab("Time") + ylab("Height (m)") +
  scale_x_datetime(breaks = seq(
    from = min(as.POSIXct(t)),
    to = max(as.POSIXct(t)),
    by = 3 * 24 * 60 * 60
  ),
  date_labels = "%d-%m") +
  theme(text = element_text(family ="serif", face="bold", size = 14)) +
  scale_color_manual(name = "Height series", values = c("GNSS-based" = "darkblue", "Atmospheric" = "red"))

###############################################
###############################################
ggplot() +
  geom_line(size = 2, aes(y = a$y, x = as.POSIXct(t), colour = "GNSS-based")) +
  geom_point(size = 2, aes(
    y = sensor2$height,
    x = as.POSIXct(sensor2$time),
    colour = "Sensor"
  )) +
  xlab("Time") + ylab("Height (m)") +
  scale_x_datetime(breaks = seq(
    from = min(as.POSIXct(t)),
    to = max(as.POSIXct(t)),
    by = 3 * 24 * 60 * 60
  ),
  date_labels = "%d-%m") +
  theme(text = element_text(family ="serif", face="bold", size = 14)) +
  scale_color_manual(name = "Height series", values = c("GNSS-based" = "darkblue", "Sensor" = "red"))




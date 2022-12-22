# sensor2
# res$height = res$height - 0.3
# res = filter (res, (day(time) >= 15 & month(time) ==9) | ( month(time) == 10))
# sensor2 = filter (sensor2, (day(time) >= 15 & month(time) ==9) | ( month(time) == 10))
sensor2$height = sensor2$height +0.1
res = filter (res, (day(time) >= 15 & month(time) ==9) | ( month(time) == 10))

resG = filter(res, group == "G")
resR = filter(res, group == "R")
resE = filter(res, group == "E")
resC = filter(res, group == "C")

Gnssdata = res[, c("height", "time", "group")] %>% group_by(group = round(res$time, "hour")) %>%
  summarise_at(vars(height), list(h_mean_hour = mean))
Gdata = resG[, c("height", "time", "group")] %>% group_by(group = round(resG$time, "hour")) %>%
  summarise_at(vars(height), list(h_mean_hour = mean))
Rdata = resR[, c("height", "time", "group")] %>% group_by(group = round(resR$time, "hour")) %>%
  summarise_at(vars(height), list(h_mean_hour = mean))
Edata = resE[, c("height", "time", "group")] %>% group_by(group = round(resE$time, "hour")) %>%
  summarise_at(vars(height), list(h_mean_hour = mean))
Cdata = resC[, c("height", "time", "group")] %>% group_by(group = round(resC$time, "hour")) %>%
  summarise_at(vars(height), list(h_mean_hour = mean))

Gnssdata$time = as.POSIXct(Gnssdata$group)
Gdata$time = as.POSIXct(Gdata$group)
Rdata$time = as.POSIXct(Rdata$group)
Edata$time = as.POSIXct(Edata$group)
Cdata$time = as.POSIXct(Cdata$group)

Gnssdata = inner_join(Gnssdata, sensor2, by="time")
Gdata = inner_join(Gdata, sensor2, by="time")
Rdata = inner_join(Rdata, sensor2, by="time")
Edata = inner_join(Edata, sensor2, by="time")
Cdata = inner_join(Cdata, sensor2, by="time")


Gnssdata = filter(Gnssdata, abs(h_mean_hour - height) < 0.25)
Gdata = filter(Gdata, abs(h_mean_hour - height) < 0.25)
Rdata = filter(Rdata, abs(h_mean_hour - height) < 0.25)
Edata = filter(Edata, abs(h_mean_hour - height) < 0.25)
Cdata = filter(Cdata, abs(h_mean_hour - height) < 0.25)

require(Metrics)
nrow(Gnssdata)
nrow(Gdata)
nrow(Rdata)
nrow(Edata)
nrow(Cdata)
cor(Gnssdata$h_mean_hour, Gnssdata$height)
cor(Gdata$h_mean_hour, Gdata$height)
cor(Rdata$h_mean_hour, Rdata$height)
cor(Edata$h_mean_hour, Edata$height)
cor(Cdata$h_mean_hour, Cdata$height)
rmse(Gnssdata$h_mean_hour, Gnssdata$height)
rmse(Gdata$h_mean_hour, Gdata$height)
rmse(Rdata$h_mean_hour, Rdata$height)
rmse(Edata$h_mean_hour, Edata$height)
rmse(Cdata$h_mean_hour, Cdata$height)


plot(Gnssdata$time, Gnssdata$h_mean_hour, type = "l", col = "darkred")
lines(Gnssdata$time, Gnssdata$height, col = "darkblue")
plot(Gdata$time, Gdata$h_mean_hour, type = "l", col = "darkred")
lines(Gdata$time, Gdata$height, col = "darkblue")


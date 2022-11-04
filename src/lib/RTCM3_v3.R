getbitu <- function(data, pos, len) {
  bit = 0.0
  for (i in (pos - 1):(pos + len - 2)) {
    bit = bitwShiftL(bit, 1) + ifelse (bitwAnd(as.numeric(data[i %/% 8 + 1]), bitwShiftR(128, i %%
                                                                                           8)) > 0x0, 1, 0)
  }
  return(bit)
}
getbitl <- function(data, pos, len) {
  bit = 0
  t = 1
  for (i in (pos - 1):(pos + len - 2)) {
    bit[t] = ifelse (bitwAnd(as.numeric(data[i %/% 8 + 1]), bitwShiftR(128, i %%
                                                                         8)) > 0x0, 1, 0)
    t = t + 1
  }
  return(bit)
}
parse.msm4 <- function(data) {
  i = 25
  type = getbitu(data, i, 12)
  if (!(type %in% c(1074, 1084, 1094, 1114, 1124))) return(NA)
  i = i + 12
  stid = getbitu(data, i, 12)
  i = i + 12
  epoch.time = getbitu(data, i, 30)
  i = i + 30
  i = i + 1 + 3 + 7 + 2 + 2 + 1 + 3
  sat.mask = which(getbitl(data, i, 64) == 1)
  i = i + 64
  sig.mask = which(getbitl(data, i, 32) == 1)
  i = i + 32
  cell.mask = getbitl(data, i, length(sat.mask) * length(sig.mask))
  i = i + length(cell.mask)
  
  sat.rr = NA
  sat.rr.mod = NA
  for (x in 1:length(sat.mask)) {
    sat.rr[x] = getbitu(data, i, 8)
    i = i + 8
  }
  for (x in 1:length(sat.mask)) {
    sat.rr.mod[x] = getbitu(data, i, 10)
    i = i + 10
  }
  
  fine.pseudo.range = NA
  fine.carrier.phase = NA
  cnr = NA
  for (x in 1:sum(cell.mask)) {
    fine.pseudo.range[x] = getbitu(data, i, 15)
    i = i + 15
  }
  for (x in 1:sum(cell.mask)) {
    fine.carrier.phase[x] = getbitu(data, i, 22)
    i = i + 22
  }
  i = i + 5 * sum(cell.mask)
  for (x in 1:sum(cell.mask)) {
    cnr[x] = getbitu(data, i, 6)
    i = i + 6
  }
  
  
  prn = rep(sat.mask, each = length(sig.mask))
  sig = rep(sig.mask, length(sat.mask))
  
  pseudoRange = NA
  phaseRange = NA
  snr = NA
  j = 1
  cLight = 299792458
  for (i in 1:length(cell.mask)) {
    if (cell.mask[i] == 1) {
      pseudoRange[i] = (cLight / 1000) * (sat.rr[(i - 1) / length(sig.mask) + 1] +
                                            sat.rr.mod[(i - 1) / length(sig.mask) + 1] / 1024 +
                                            2 ^ (-24) * fine.pseudo.range[j])
      phaseRange[i] = (cLight / 1000) * (sat.rr[(i - 1) / length(sig.mask) + 1] +
                                           sat.rr.mod[(i - 1) / length(sig.mask) + 1] / 1024 +
                                           2 ^ (-29) * fine.carrier.phase[j])
      snr[i] = cnr[j]
      j = j + 1
    }
    else {
      pseudoRange[i] = NA
      phaseRange[i] = NA
      snr[i] = NA
    } 
  }
  
  if(!dir.exists(".temp/rtcm/")) {
    dir.create(".temp/rtcm/")
  }
  if (!file.exists(paste0(".temp/rtcm/tempobs.", type))) {
    write.table(tidyr::drop_na(data.frame(prn, sig, pseudoRange, phaseRange, snr, epoch.time)), paste0(".temp/rtcm/tempobs.", type), row.names = F, col.names = T, append =T)
  }
  else{
    write.table(tidyr::drop_na(data.frame(prn, sig, pseudoRange, phaseRange, snr, epoch.time)), paste0(".temp/rtcm/tempobs.", type), row.names = F, col.names = F, append = T)
  }
}

cut.out <- function(data) {
  i = 2
  pb <- progress_bar$new(
    format = "  Parsing RTCM [:bar] :percent ETA: :eta eslapsed: :elapsed",
    total = length(data))
  while (i < length(data)) {
    len = bitwShiftL(as.numeric(data[i]), 8) + as.numeric(data[i + 1])
    parse.msm4(data[(i - 1):(i + len + 3)])
    i <- i + len + 6
    pb$tick(i-(i-len-6))
  }
  gc()
}

file = "input/2143-4.dat"
bin = readBin(file, "raw" , file.info(file)$size)
dat = cut.out(bin)

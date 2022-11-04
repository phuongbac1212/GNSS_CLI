#include<Rcpp.h>
using namespace Rcpp;


unsigned int getbitu(const unsigned char *buff, int pos, int len)
{
  unsigned int bits=0;
  int i;
  for (i=pos;i<pos+len;i++) bits=(bits<<1)+((buff[i/8]>>(7-i%8))&1u);
  return bits;
}

DataFrame parse_msm4(const unsigned char *data) {
  int i = 25;
  int type = getbitu(data, i, 12);
  i += 12;
  int stid = getbitu(data, i, 12);
  i = i + 12;
  int epoch = getbitu(data, i, 30);
  i = i + 30;
  i = i + 1 + 3 + 7 + 2 + 2 + 1 + 3;
  int tmp = getbitu(data, i, 64);
  NumericVector satMask;
  for (int i =1; i<=64; i++) {
    if (tmp >> (64-i) == 1) {
      satMask.push_back(i);
    }
  }
  
  i = i + 64;
  tmp = getbitu(data, i, 32);
  NumericVector sigMask;
  
    for (int i =1; i<=32; i++) {
      if (tmp >> (32-i) == 1) {
        sigMask.push_back(i);
      }
    }
  i = i + 32;
  tmp = getbitu(data, i, satMask.length() * sigMask.length());
  NumericVector cellMask;
  for (int i =1; i<= satMask.length() * sigMask.length(); i++) {
    if (tmp >> ( satMask.length() * sigMask.length()-i) == 1) {
      cellMask.push_back(i);
    }
  }
  i = i + cellMask.length();
  /*
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
   */
      return(DataFrame::create(type, epoch));
}
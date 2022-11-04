#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <cstdio>
using namespace Rcpp;

// [[Rcpp::export]]
bool skytraq_checksum(IntegerVector msg) {
  int len = msg.length(), i= 4, cs = 0;
  for (; i<len-3; i++) {
    cs ^= int(msg[i]);
  }
  return (cs == msg[len - 3]);
}

// // [[Rcpp::export]]
// String skytraq_getObsType(int type, int svrID) {
//   String name;
//   switch (type) {
//   case 0:
//     name = "PG"; break;
//   case 1:
//     name = "PS"; break;
//   case 2: 
//     name= "PR";break;
//   case 3:
//     name = "PE"; break;
//   case 4:
//     name = "PJ"; break;
//   case 5:
//     name = "PC";break;
//   case 6:
//     name = "PI"; break;
//   }
//   
//   int id;
//   switch (type) {
//   case 0:
//     id = svrID; break;
//   case 1:
//     id = svrID -120 ; break;
//   case 2: 
//     id = svrID ; break;
//   case 3:
//     id = svrID ; break;
//   case 4:
//     id = svrID ; break;
//   case 5:
//     id = svrID -192 ; break;
//   case 6:
//     id = svrID ; break;
//   }
//   char obs[10];
//   sprintf(obs, "%s%2d", name.get_cstring(), id);
//   return (toString(obs));
// }

// [[Rcpp::export]]
List skytraq_parseE5(IntegerVector msg , String path) {
  Function determine_lambda("determine_lambda");
  Function skytraq_getObsType("skytraq.getObsType");
  int i = 4;
  int MsgID = msg[i++];
  int Version = msg[i++];
  int IOD = msg[i++];
  int WN = msg[i]<<8|msg[i+1]; i+=2;
  int TOW = msg[i] << 24 | msg[i+1] << 16 | msg[i+2] <<8 | msg[i+3]; i+=4;
  int Measurement_period = msg[i]<<8|msg[i+1]; i+=2;
  int Measurement_indicator = msg[i]; i+=2;
  int NMEAS = msg[i];
  List l = List::create();
  
  int GNSS_TYPE, Signal_type, SVID, FreqIDnLTI, CNR, Pseudorange_SD, carrier_cycle_SD, doppler_SD, Channel_Indicator;
  long long Pseudorange, carrier_phase, Doppler;
  
  char outfile[200];
  sprintf(outfile, "%s/%d.out", path.get_cstring(), WN * 7 * 24 * 60 * 60 + TOW / 1000);
  
  std::ofstream fout;
  fout.open(outfile, std::ios::app);
  
  fout << "GNSS Type, Signal Type, Server ID, Wavelength, CNR, Pseudorange, Accumulated carrier cycle, Doppler frequency, Observation Type, Timestamp" << std::endl;
  for (int loop = 0; loop<NMEAS; loop++){
    i++;
    GNSS_TYPE = msg[i] & 7;
    Signal_type = (msg[i] >> 4) & 0xF; i++;
    SVID = msg[i++];
    FreqIDnLTI = msg[i++];
    CNR = msg[i++];
    Pseudorange = msg[i++];
    for (int j =7; j>0; j--) {Pseudorange = Pseudorange << 8 | msg[i++];}
    carrier_phase = msg[i++];
    for (int j =7; j>0; j--) {carrier_phase = carrier_phase << 8 | msg[i++];}
    Doppler = msg[i++];
    for (int j =3; j>0; j--) {Doppler = Doppler << 8 | msg[i++];}
    Pseudorange_SD = msg[i++];
    carrier_cycle_SD = msg[i++];
    doppler_SD = msg[i++];
    Channel_Indicator = msg[i];
    i += 3;
      
    fout << std::fixed <<
      GNSS_TYPE << "," <<
      Signal_type << "," <<
      SVID << "," <<
      *REAL(determine_lambda(GNSS_TYPE, Signal_type, (FreqIDnLTI >> 4) - 7) )<< "," <<
      CNR << "," <<
      *(double*) &Pseudorange << ","  <<
      *(double*) &carrier_phase << "," <<
      *(float*) &Doppler << "," <<
      Rcpp::as<std::string>(skytraq_getObsType(GNSS_TYPE, SVID)) << "," <<
      WN * 7 * 24 * 60 * 60 + TOW / 1000 << std::endl;
    
    // DataFrame df = DataFrame::create(
    //   Named("GNSS Type") = GNSS_TYPE,
    //   Named("Signal Type") = Signal_type,
    //   Named("Server ID") = SVID,
    //   Named("Wavelength") = determine_lambda(GNSS_TYPE, Signal_type, (FreqIDnLTI >> 4) - 7),
    //   Named("CNR") = CNR,
    //   Named("Pseudorange") = *(double*) &Pseudorange,
    //   Named("Accumulated carrier cycle") = *(double*) &carrier_phase,
    //   Named("Doppler frequency") = *(float*) &Doppler,
    //   Named("Observation Type") = skytraq_getObsType(GNSS_TYPE, SVID),
    //   Named("Timestamp") = WN * 7 * 24 * 60 * 60 + TOW / 1000
    // );
    // l.push_back(df);
  }
  i++;
  return(l);
}

/*** R
# skytraq.parseE5 <- function(msg, path = parser_path) {
#    s = skytraq_parseE5(msg) %>% data.table::rbindlist()
#   write_fst(s, paste0(path, unique(s$Timestamp), ".fst"))
# }
# s = skytraq_parseE5(msg, "/home/fang/Working/GNSS-R CLI/output/parser") %>% data.table::rbindlist()
*/


### some definite path
input_path = "input/"
input_backup_path = paste0(input_path, "/backup/")
extract_path = "output/extract/"
extract_backup_path  = "output/extract/backup/"
parser_path = "output/parser/"
result_path = "output/result/"
result_backup_path = "output/result/backup/"
shp_path = "input/shp/"
sp3_path = "input/sp3/"
sp3_ftp = "ftp://igs.ign.fr/pub/igs/products/mgex/"
#sp3_ftp = "ftp://gdc.cddis.eosdis.nasa.gov/pub/gps/products/mgex/"

### some global value for extract
sp3_suffix = "WUM0MGXULA_"
gps_time_offset = 315964782
speed_of_light = 299792458
n.cluster = 6
rec = c(x = -1845476.391,
        y = 5828929.663,
        z = 1810098.982)

elevation_min = 2 * pi/180
elevation_max = 25 * pi/180
# elevation_min = 30 * pi/180
# elevation_max = 70 * pi/180
azimuth_to_take = c(150, 250)* pi/180

## some global value for analysis
data_series_length = 1200
minimum_length_series = 500
window.step = 600
max_hole = 20
plot_snr = F
plot_lsp = F
h_fixed = FALSE

h0 = 7.84
h_max = 9.696
h_min = 5.985
h_point_max = 1e-5

method = "Larson"

### some
do_extract = F
do_analysis = T
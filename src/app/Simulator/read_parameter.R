##################################################################
#
# janvier 2016 - Nicolas ROUSSEL
#
##################################################################
library(data.table)
library(shinyFiles)
library(shiny)
param = reactiveValues()
#read parameter before all
volumes <-
  c(Home = fs::path_home(),
    "R Installation" = R.home(),
    getVolumes()())
parameter.read.raw = function() {
  file = "cfg/parameters_temp.tmp"
  param <<-
    scan(
      file,
      sep = ':',
      na.strings = '-',
      what = character(0),
      quiet = TRUE,
      blank.lines.skip = TRUE,
      comment.char = '#'
    )
  param <<- matrix(param, ncol = 2, byrow = TRUE)
  params = as.data.frame(param)
  param = params$V2
  names(param) = params$V1
  save(param, file = "cfg/parameters_temp.RDS")
  return(param)
}

parameter.update.from.input <- function(input) {
  try(expr = {param["dem_file_path"] <<- parseFilePaths(roots = c(wd = '~'), input$DEM_file)$datapath}, silent=T)
  try(expr = {param["algorithm"] <<- input$PARAM_algo}, silent=T)
  try(expr = {param["basemap"] <<- input$PARAM_basemap}, silent=T)
  
  try(expr = {param["low_horizontal_mask"] <<- input$PARAM_visibility_hor[1]}, silent=T)
  try(expr = {param["up_horizontal_mask"] <<- input$PARAM_visibility_hor[2]}, silent=T)
  try(expr = {param["low_vertical_mask"] <<- input$PARAM_visibility_ver[1]}, silent=T)
  try(expr = {param["up_vertical_mask"] <<- input$PARAM_visibility_ver[2]}, silent=T)
  try(expr = {param["wavelength"] <<- input$PARAM_wavelength}, silent=T)
  try(expr = {param["receiver_eheight"] <<- input$REC_eh}, silent=T)
  try(expr = {param["receiver_height"] <<- input$REC_height}, silent=T)
  try(expr = {param["receiver_latitude"] <<- input$input$REC_lat}, silent=T)
  try(expr = {param["receiver_longitude"] <<- input$REC_long}, silent=T)
  try(expr = {param["sp3_file_path"] <<- parseFilePaths(roots =volumes ,input$SAT_file)$datapath}, silent=T)
  try(expr = {param["SP3_limit_begin"] <<- input$SP3_limit_begin}, silent=T)
  try(expr = {param["SP3_limit_end"] <<- input$SP3_limit_end}, silent=T)
  try(expr = {param["SP3_limit_end"] <<- input$SP3_limit_end}, silent=T)
  try(expr = {param["constelations"] <<- noquote(paste(input$constelations, collapse = ','))}, silent=T)
    # print(parseFilePaths(roots =volumes ,input$SAT_file))
  save(param, file = "cfg/parameters_temp.RDS")
}

dt_output = function(title, id) {
  fluidRow(column(12, h1(paste0(
    'Table ', sub('.*?([0-9]+)$', '\\1', id), ': ', title
  )),
  hr(), DTOutput(id)))
}
render_dt = function(data,
                     editable = 'cell',
                     server = TRUE,
                     ...) {
  renderDT(
    data,
    selection = 'none',
    server = server,
    editable = editable,
    ...
  )
}


render_param_dt = function(params,
                           editable = 'cell',
                           server = TRUE,
                           ...) {
  parameters = data.table(names(params), params)
  renderDT(
    parameters,
    selection = 'none',
    server = server,
    editable = editable,
    ...
  )
}

tryCatch(
  load("cfg/parameters_temp.RDS"),
  error = function(con) {
    print(paste("Cannot open parameter file [", con, "]"))
    param = parameter.read.raw()
    print(paste("reading [", length(param), "] variables from raw file"))
  }
)

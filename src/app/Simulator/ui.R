#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#setwd("src/app/Simulator/")

library(shiny)
require(shinyFiles)
require(shinyWidgets)
require(leaflet)
library(DT)



source("read_parameter.R")
#read parameter before all

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  #includeCSS("assert/css/dark_mode.css"),
  shinyjs::useShinyjs(),
  htmltools::htmlDependency(
    name = "font-awesome",
    version = "99.0",
    src = "assert/fontawesome",
    stylesheet = "css/all.min.css"
  ),
  # Application title
  titlePanel(h1(strong("GNSS-R Simulator"))),
  
  # Sidebar with a slider input for number of bins
  tabsetPanel(
    id = "tabset_menu",
    tabPanel(
      "Settings",
      h1(strong("Digital Elevation Model")),
      checkboxInput(inputId = "checkbox_DEM", label = "DEM intergrate or not ?"),
      uiOutput("ui_DEM"),
      h1(strong("Receiver Position")),
      numericInput("REC_height", "Receiver height(m)", param["receiver_height"]),
      h4("Receiver"),
      fluidRow(
        column(4, numericInput("REC_long", 'longitude', param["receiver_longitude"])),
        column(4, numericInput("REC_lat", "latitude", param["receiver_latitude"])),
        column(4, numericInput("REC_eh", "elevation", param["receiver_eheight"]))
      ),
      br(),
      fluidRow(
        actionButton("SP3_down", "Down a new SP3"),
        shinyFilesButton(
          'SAT_file',
          'Select ephemeris file',
          'Please select a .SP3, .alm or .pos file',
          FALSE,
          icon = icon("satellite")
        ),
      ),
      h1(strong("Program parameters")),
      h4("Satellite visibility"),
      numericRangeInput(
        "PARAM_visibility_ver",
        "Vertical Field of visibility",
        value = c(param["mask_ver_min"], param["mask_ver_max"]),
        min =0, max = 90
      ),
      numericRangeInput(
        "PARAM_visibility_hor",
        "Horizontal Field of visibility",
        value = c(param["mask_hor_min"], param["mask_hor_max"]), 
        min=0, max = 360
      ),
      radioButtons(
        "PARAM_algo",
        "Algorithm used for calculation",
        selected = param["algorithm"],
        choices = c("ellipsoid",
                    "sphere",
                    "plane")
      ),
      icon = icon("gear", lib = "font-awesome"),
      h4("Satellite visibility"),
      radioButtons(
        "PARAM_wavelength",
        "Wavelength",
        selected = param["wavelength"],
        choices = c("L1", "L2", "L5")
      ),
      actionButton("BTN_start", "Simulate", icon = icon("robot")),
      br(),
      print("~~~")
      ,
      value = "tab1"
    ),
    tabPanel(
      "Reflected Point",
      leafletOutput("reflect_point_map"),
      value = "tab2",
      icon = icon("explosion")
    ),
    tabPanel(
      "Feshnel Zone",
      leafletOutput("feshnel_zone"),
      value = "tab3",
      icon = icon("earth-asia")
    ),
    tabPanel(
      "Parameters",
      dt_output("Parameter for simulator",
                'paramter_tab'),
      fluidRow(
        actionButton("BTN_param_ok", "OK", icon = icon("check")),
        actionButton("BTN_param_close", "Close", icon = icon("xmark"))
      )
      ,
      value = "tab4"
    )
  )
))

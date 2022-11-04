library(shiny)
require(shinyFiles)
require(leaflet)
require(raster)

source("read_parameter.R")

shinyServer(function(input, output, session) {
  # save(input, file= "cfg/input.RDS")
  
  output$ui_DEM <- renderUI(if (input$checkbox_DEM) {
    fluidPage(
      shinyFilesButton(
        'DEM_file',
        'Select DEM file',
        'Please select a DEM file',
        FALSE,
        icon = icon("mountain-city")
      ),
      p(
        strong("NOTE"),
        em(
          "The greater part of the DEM is, the longer the calculation would be. Recommend: +/- 2 deg"
        )
      ),
      numericInput("DEM_sample_interval", 'Sampling interval:', 1),
      p("x 0.00083 deg"),
      radioButtons(
        inputId = "DEM_iterpoolation_type",
        label = "Iterpolation:",
        choices = c("cubic", "bilinear")
      ),
      fluidRow(column(
        6, actionButton("DEM_export_kml", "Create a KML file")
      ), column(
        6, actionButton("DEM_display", "Display DEM")
      ))
    )
    #numericInput()
  })
  
  shinyFileChoose(
    input,
    'DEM_file',
    roots = c(wd = '~'),
    #filetypes = c('', 'txt'),
    defaultPath = '',
    defaultRoot = 'wd'
  )
  
  shinyFileChoose(
    input,
    'SAT_file',
    roots = c(wd = '~'),
    filetypes = c('sp3', 'SP3', "alm", "pos"),
  )
  
  
  
  observeEvent(input$DEM_display, {
    showModal(modalDialog(
      title =  "Visualizing DEM",
      footer = modalButton("Close"),
      renderPlot({
        if (input$checkbox_DEM) {
          file = parseFilePaths(roots = c(wd = '~'), input$DEM_file)
          tryCatch(
            expr = {
              dem.tif = raster(file$datapath)
            },
            error = function(con) {
              showNotification(paste("Cannot open DEM file [", con, "]"))
            }
          )
          
        }
        else {
          showNotification("Please select DEM file")
        }
        plot(dem.tif)
      })
    ))
  })
  
  # function to handle start buttom
  reflect_out <- eventReactive(input$BTN_start, {
    #updateTabsetPanel(session, "tabset_menu", 'tab2')
    source("shiny_principal.R")
    require(randomcoloR)
    pos_S_geo = principal(input)
    colorPalette = randomColor(length(unique(pos_S_geo$ID)))
    names(colorPalette) = unique(pos_S_geo$ID)
    pos_S_geo$color = unlist(lapply(pos_S_geo$ID, FUN = function(x) {colorPalette[x]}))
    leaflet() %>% 
      addProviderTiles(provider = providers$Esri.WorldImagery) %>% 
      addCircleMarkers(lng = pos_S_geo[,1], lat = pos_S_geo[,2], radius = 1, color = pos_S_geo$color) %>% 
      setView(lng = as.numeric(param["receiver_longitude"]), lat = as.numeric(param["receiver_latitude"]), zoom = 16)
  })
  
  output$reflect_point_map <- renderLeaflet(reflect_out())
  
  output$paramter_tab = render_param_dt(param, "cell", options = list(pageLength = 100))
  
  observeEvent(input$BTN_start, {
    updateTabsetPanel(session, "tabset_menu", 'tab2')
    
  })
  observeEvent(input$paramter_tab_cell_edit, {
    params = data.table(names(param), param)
    
    info = input$paramter_tab_cell_edit
    params <- editData(params, info)
    temp = c(params$param)
    names(temp) = params$V1
    param <<- temp
    
    showNotification("Parameters cell edited")
  })
  observeEvent(input$BTN_param_ok, {
    save(param, file = "cfg/parameters_temp.RDS")
    showNotification("Parameters updated")
    updateTabsetPanel(session, "tabset_menu", 'tab1')
  })
  observeEvent(input$BTN_param_close, {
    updateTabsetPanel(session, "tabset_menu", 'tab1')
  })
  observeEvent(input$tabset_menu, {
    if (input$tabset_menu == "tab4") {
      parameter.update.from.input(input)
      save(param, file = "cfg/parameters_temp.RDS")
    }
  })
})

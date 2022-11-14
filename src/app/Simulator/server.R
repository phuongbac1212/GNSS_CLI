library(shiny)
require(shinyFiles)
require(leaflet)
require(raster)

source("read_parameter.R")
volumes <-
  c(Home = fs::path_home(),
    "R Installation" = R.home(),
    getVolumes()())
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
    defaultPath = '~',
    defaultRoot = 'wd'
  )
  
  shinyDirChoose(input,
                 "SP3_dir_btn",
                 roots = volumes,
                 session = session)
  
  
  shinyFileChoose(
    input,
    "SAT_file",
    roots = volumes,
    session = session,
    filetypes = c('sp3', 'SP3', "alm", "pos")
  )
  
  observeEvent(input$SP3_down, {
    showModal(
      modalDialog(
        title =  "Donwload SP3",
        footer = tagList(
          modalButton("Cancel"),
          actionButton("SP3_start_down", "OK")
        ),
        dateInput("SP3_date", label = "Choose a date to down!"),
        shinyDirButton(
          "SP3_dir_btn",
          "Choose a place to save",
          "Select where to save!"
        ),
        textOutput("SP3_console")
      )
    )
  })
  
  
  
  observeEvent(input$SP3_start_down, {
    withCallingHandlers({
      shinyjs::html(id = "SP3_console", html = "")
      source("sp3.R")
      time = as.numeric(as.POSIXct(input$SP3_date))
      sp3_dir = parseDirPath(volumes , input$SP3_dir_btn)
      
      print(time)
      print(sp3_dir)
      
      link = sp3.getLink(timestamp = 1668038400 - as.numeric(param["gps_time_offset"]))
      sp3.download(link = link, path = sp3_dir)
    },
    message = function(m) {
      shinyjs::html(id = "SP3_console",
                    html = m$message,
                    add = TRUE)
    },
    warning = function(m) {
      shinyjs::html(id = "SP3_console",
                    html = m$message,
                    add = TRUE)
    })
    #
    # output$SP3_console <- renderPrint({
    #   source("sp3.R")
    #   time = as.numeric(as.POSIXct(input$SP3_date))
    #   sp3_dir = parseDirPath(volumes , input$SP3_dir_btn)
    #
    #   print(time)
    #   print(sp3_dir)
    #
    #   link = sp3.getLink(timestamp = 1668038400- as.numeric(param["gps_time_offset"]))
    #   sp3.download(link = link, path = sp3_dir)
    #   # sp3.download(sp3.getLink(as.numeric(input$SP3_date))
    # })
  })
  
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
    res = principal(input)
    return(res)
  })
  
  output$reflect_point_map <- renderLeaflet({
    res = reflect_out()
    map = leaflet() %>% addTiles()
    map %>% addCircleMarkers(lng = res$lon, lat = res$lat, radius = 1)
  })
  output$feshnel_zone <- renderLeaflet({
    res = reflect_out()
    map = leaflet() %>% addTiles()
    ellipse = apply(res, 1, function(p, map) {
      # print(p)
      R = 637100
      
      dLat = as.numeric(p["a"]) / R
      dLon = as.numeric((p["b"])) / (R * cos(pi * as.numeric(p["lat"]) /
                                               180))
      ell = DrawEllipse(
        x = as.numeric(p["lon"]),
        y = as.numeric(p["lat"]),
        radius.y = dLon,
        radius.x = dLat,
        rot = as.numeric(p["azimuth"]),
        plot = F
      )
      return(list(
        lng = ell$x,
        lat = ell$y,
        ID = p["ID"]
      ))
      
    })
    
    for (e in ellipse) {
      e$lng[length(e$lng) + 1] = e$lng[1]
      e$lat[length(e$lat) + 1] = e$lat[1]
      map = map %>% addPolygons(lng = e$lng, lat = e$lat)
    }
    map
  })
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

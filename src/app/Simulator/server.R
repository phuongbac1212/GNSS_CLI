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
    pos_S_xyz = res[, c("XS", "YS", "ZS")]
    pos_S_geo = data.frame(t(apply(pos_S_xyz, 1, pos2geo, R_0 = as.numeric(param["semi_major_WGS84"]), R_p = as.numeric(param["semi_minor_WGS84"]))))
    names(pos_S_geo) = c("long", "lat", "height")
    map = leaflet() %>% addProviderTiles(providers[[param[["basemap"]]]])

    map %>% addCircleMarkers(
      lng = pos_S_geo$long,
      lat = pos_S_geo$lat,
      radius = 1,
      color = res$color
    )

  })
  output$feshnel_zone <- renderLeaflet({
    res = reflect_out()
    zone = apply(res, 1, function(x) {
      a <- x[["a"]]
      b <- x[["b"]]
      aire <- x[["c"]]
      
      nb_gones <-
        40 # How many sides of the ellipse (= approximated as a polygone)
      
      nb_gones <- nb_gones + 2
      # Coordinates of the known points in the 3D WGS84 system
      
      x1 <- x[["XA"]]
      y1 <- x[["YA"]]
      z1 <- x[["ZA"]]
      
      x2 <- x[["XB"]]
      y2 <- x[["YB"]]
      z2 <- x[["ZB"]]
      
      x3 <- x[["XAABB"]]
      y3 <- x[["YAABB"]]
      z3 <- x[["ZAABB"]]
      
      x0 <-  x[["XS"]]
      y0 <-  x[["YS"]]
      z0 <-  x[["ZS"]]
      
      # Coordinates of the known points in the 3D local system, defined by X : axis S->A, Y : axis S->B, Z : the straight line orthogonal to forme a direct trihedral.
      
      X1 <- a
      Y1 <- 0
      Z1 <- 0
      
      X2 <- 0
      Y2 <- b
      Z2 <- 0
      
      X3 <- -a
      Y3 <- 0
      Z3 <- -b
      
      X0 <- 0
      Y0 <- 0
      Z0 <- 0
      
      u <- 0
      v <- 0
      
      source("geo2pos.R")
      pos_S_xyz <- as.numeric(c(x0, y0, z0))
        # geo2pos(R_0_WGS84, R_p_WGS84, c(mon_data[num_ligne, 1], mon_data[num_ligne, 2], mon_data[num_ligne, 3]))
      
      # We look for the parameters to go from a system to the other.
      # To do that, we do a 9 parameters transformation :
      # We have :
      # x = T1 + A*X + B*Y + C*Z
      # y = T2 + D*X + E*Y + F*Z
      # z = T3 + G*X + H*Y + I*Z
      # and so we have: A %*% X = B, with :
      # A =   X1 Y1 Z1 0 0 0 0 0 0
      #       0 0 0 X1 Y1 Z1 0 0 0
      #       0 0 0 0 0 0 X1 Y1 Z1
      #       X2 Y2 Z2 0 0 0 0 0 0
      #       0 0 0 X2 Y2 Z2 0 0 0
      #       0 0 0 0 0 0 X2 Y2 Z2
      #       X3 Y3 Z3 0 0 0 0 0 0
      #       0 0 0 X3 Y3 Z3 0 0 0
      #       0 0 0 0 0 0 X3 Y3 Z3
      # X =  c(A,B,C,D,E,F,G,H,I)
      # B = (x1,y1,z1,x2,y2,z2,x3,y3,z3)
      
      t1 <- pos_S_xyz[1]
      t2 <- pos_S_xyz[2]
      t3 <- pos_S_xyz[3]
      
      A <-
        matrix(
          c(
            X1,
            Y1,
            Z1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            X1,
            Y1,
            Z1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            X1,
            Y1,
            Z1,
            X2,
            Y2,
            Z2,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            X2,
            Y2,
            Z2,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            X2,
            Y2,
            Z2,
            X3,
            Y3,
            Z3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            X3,
            Y3,
            Z3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            X3,
            Y3,
            Z3
          ),
          ncol = 9,
          byrow = TRUE
        )
      B <- c(x1 - t1,
             y1 - t2,
             z1 - t3,
             x2 - t1,
             y2 - t2,
             z2 - t3,
             x3 - t1,
             y3 - t2,
             z3 - t3)
      
      result <- qr.solve(A, B)
      A <- result[1]
      B <- result[2]
      C <- result[3]
      D <- result[4]
      E <- result[5]
      F <- result[6]
      G <- result[7]
      H  <- result[8]
      I <- result[9]
      
      # XX : abscissa of the differents points of the ellipse we will plot
      XX1 <- seq(
        from = u - a,
        to = u + a,
        length = round(nb_gones / 2)
      )
      XX2 <- seq(
        from = u + a,
        to = u - a,
        length = round(nb_gones / 2)
      )
      XX <- c(XX1, XX2)
      # YY : ordinates of the differents points of the ellipse we will plot
      YY1 <- sqrt(b ^ 2 - ((XX1 - u) ^ 2 / a ^ 2 * b ^ 2)) + v
      YY2 <- -sqrt(b ^ 2 - ((XX2 - u) ^ 2 / a ^ 2 * b ^ 2)) + v
      YY <- c(YY1, YY2)
      ZZ <- seq(from = 0,
                to = 0,
                length = length(YY))
      xx <- 0
      yy <- 0
      zz <- 0
      long <- 0
      lat <- 0
      he <- 0
      # We transform the points coordinates of the ellipse we will plot from the 3D local system to the 3D WGS84 system thanks to the 9 parameters we determined before.
      for (u in 1:length(XX)) {
        xx[u] <- t1 + A * XX[u] + B * YY[u] + C * ZZ[u]
        yy[u] <- t2 + D * XX[u] + E * YY[u] + F * ZZ[u]
        zz[u] <- t3 + G * XX[u] + H * YY[u] + I * ZZ[u]
        result <-
          pos2geo(as.numeric(param["semi_major_WGS84"]),
                  as.numeric(param["semi_minor_WGS84"]), 
                  c(xx[u], yy[u], zz[u]))
        long[u] <- result[1]
        lat[u] <- result[2]
        he[u] <- result[3]
      }
      if (x2 == y2 &
          y2 == z2 &
          z2 == 0) {
        # If an error occured during the calculation of B (sometimes, the discriminant is negative, I don't know why) --> we don't plot this Fresnel surface
        long <- x0
        lat <- y0
      }
      return(list(long, lat, x[["color"]]))
    })
    map = leaflet() %>% addProviderTiles(providers[[param[["basemap"]]]])
    for (z in zone) {
      map = map %>% addPolygons(lng = z[[1]], lat = z[[2]], fillColor = as.character(z[[3]]), fillOpacity = 1, opacity = 0.5)
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

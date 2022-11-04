#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#sink(paste("logs/", Sys.time(), ".txt", sep= ""), append=FALSE, split=TRUE)

require(shiny)
require(shinyWidgets)
require(fst)
require(DT)
# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    title = "Đọc FST",
    column(
        2,
        fileInput(inputId = "file",
                  label = "File: "),
        submitButton("Update View", icon("refresh")),
    ),
    hr(),
    dataTableOutput("table")
    
))
server <- shinyServer(function(input, output) {
    output$table <- renderDataTable({
        print(input$file)
        if (is.null(input$file)) {
            print("no data")
            return(NULL)
        }
        res = (as.data.frame(fst(input$file$datapath)))
        print(head(res))
        return(res)
    })
    
})

shinyApp(ui, server)

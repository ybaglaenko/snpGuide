library(shiny)
library(DT)
source("snpGuide_global.R")
version = "0.1"
update = "May 30th, 2024"
changes = "Going live ... "


ui <- fluidPage(
  pageWithSidebar(
    
    # App title ----
    headerPanel(paste0("snpGuide v", version)),
    # Sidebar panel for inputs ----
    sidebarPanel(
      p(""),
      em(paste0("Updated ", update, ".")),
      p(""),
      em(paste0("Lastest Change: ", changes, ".")),
      p(""),
      em(paste0("Please email yuriy.baglaenko@cchmc.org with any questions, comments, or concerns.")),
      p(""),
      tags$a(href="https://www.baglaenkolab.com", em("Check us out at baglaenkolab.com")),
      p(""),
      textInput("rsID", "Input rsID to target"),
      p(""),
      fileInput("upload", "Upload batch rsIDs as .txt file", accept = ".txt"),
        
      textInput("pam", "PAM sequence", value = "NGG"),    
      
      sliderInput(inputId = "editing_window",
                  label = "Choose base editing window",
                  min = 1,
                  max = 20,
                  value = c(3,10),
                  step = 1)
      
    ),
    
    # Main panel for displaying outputs ----
     mainPanel(
          h2("Predicted guides"),
          br(),
          DT::DTOutput("df"),
          br(),
          downloadButton("download_guide", "Download all guides"))
    ))

###### Server side function
server <- function(input, output, session) {

guides <- reactive({
  if (is.null(input$upload)) {
    # Single rsID mode
    snpGuide(
      rsid = as.character(input$rsID),
      pam = input$pam,
      guide_length = 20,
      min_editing_window = input$editing_window[1],
      max_editing_window = input$editing_window[2],
      flank = 20
    )
  } else {
    # Batch mode: read lines from uploaded file
    rsids <- readLines(input$upload$datapath)
    rsids <- unique(trimws(rsids))  # clean whitespace, remove duplicates

    # Safely apply snpGuide to each rsID and combine results
    map_df(rsids, ~ possibly(snpGuide, otherwise = NULL)(
      rsid = .x,
      pam = input$pam,
      guide_length = 20,
      min_editing_window = input$editing_window[1],
      max_editing_window = input$editing_window[2],
      flank = 20
    ))
  }
})
    
### Allow downloading
output$download_guide <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"guides", input$rsID, ".csv", sep = "_")
    },
    content = function(file) {
      write.csv(guides(), file, row.names = FALSE)
    }
  )
   
### Render outputs
     output$df = DT::renderDT({
      datatable(guides())
     })

}
shinyApp(ui, server)

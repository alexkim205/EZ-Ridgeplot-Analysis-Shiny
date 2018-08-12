library(shiny)
# library(shinycssloaders)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("ScreenSeq Visualization"),
  
  # Main plot
  tabsetPanel(
    id = "mainpanel",
    tabPanel("Screen by Perturbation",
             uiOutput("mainplotui", 
                      align = "center",
                      style = "height: 60vh; 
                        overflow-y: auto;
                        margin-top: 30px;"),
             hr(),
             fluidRow(
               # Parameters
               column(3,
                      h4("Parameters"),
                      fileInput('file', 'Choose CSV Files', accept=c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv"),
                        multiple = TRUE)
                      ),
               column(3,
                      uiOutput('pert'),
                      uiOutput('pert_type'),
                      uiOutput('rep')
                      # selectInput('pert', 'Perturbation', c(1, 2)),
                      # selectInput('pert_type', 'Perturbation type', c(1, 2, 3)),
                      # selectInput('rep', 'Replicate', c(1, 2, 3, 4))
               ),
               column(3,
                      uiOutput('day'),
                      uiOutput('readout_gene')
               ),
               column(3,
                      radioButtons('plottype', 'Type of Plot', c("Ridge" = "ridge", 
                                                                 "Violin" = "violin",
                                                                 "Box" = "box")),
                      hr(),
                      actionButton('submit', "Generate plot"),
                      downloadButton('downloadPlot')
                      # bookmarkButton()
               )
             )
    ),
    tabPanel("Screen by Readout Gene",
             fluidRow(
               column(4, style = "padding: 20px;",
                      h4("Parameters"),
                      uiOutput('day_drug'),
                      uiOutput('readout_gene_drug'),
                      actionButton('submit_drug', "Generate plot"),
                      downloadButton('downloadDrugPlot')
               ),
               column(8,
                      uiOutput("drugplotui", 
                               align = "center",
                               style = "height: 150vh; 
                               overflow-y: auto;
                               margin-top: 30px;"))
             )
    )
    
  )
))

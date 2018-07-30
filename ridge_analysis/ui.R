library(shiny)
# library(shinycssloaders)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("ScreenSeq Ridge Plots"),
  
  # Main plot
  tabsetPanel(id = "mainpanel",
              tabPanel("tabPlot", align="center",
                       style = "height: 60vh; 
                       overflow-y: auto;
                       margin-top: 30px;",
                       uiOutput("mainplotui")
                       # withSpinner(
                       #   uiOutput("mainplotui")
                       # )
              )
  ),
  
  hr(),
  
  fluidRow(
    # Parameters
    column(3,
           h4("Parameters"),
           fileInput('file', 'Choose a CSV File', accept=c(
             "text/csv",
             "text/comma-separated-values,text/plain",
             ".csv"))),
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
           uiOutput('readout_gene'),
           actionButton('submit', "Generate plot"),
           downloadButton('downloadPlot'),
           bookmarkButton()
           # actionButton("submit", "Generate plot", icon("refresh"))
           # selectInput('day', 'Day', c(7, 12, 15)),
           # selectInput('gene', 'Gene', c(1, 2, 3, 4, 5, 6))
    ),
    column(3,
           radioButtons('plottype', 'Type of Plot', c("Ridge Plot" = "ridge", 
                                                      "Violin" = "violin",
                                                      "Box" = "box"))
    )
  )
))

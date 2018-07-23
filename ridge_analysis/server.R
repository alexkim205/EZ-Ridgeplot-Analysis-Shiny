library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(tidyverse)
library(ggridges)

# ── Conflicts ──
# ✖ dplyr::filter() masks stats::filter()
# ✖ dplyr::lag()    masks stats::lag()
filter <- dplyr::filter
lag <- dplyr::lag

source("plot.R") # Import make_summary_plot() function

shinyServer(function(input, output, session) {
  
  width_per_day = 250
  height_per_replicant = 20
  
  # Store reactive values
  v <- reactiveValues(
    data = NULL,
    makeplot = FALSE
  )
  
  tab_width = eventReactive(input$submit, {
    width_per_day * length(input$"_day")
  })
  tab_height = eventReactive(input$submit, {
    row_height <- length(input$"_rep") * height_per_replicant
    row_height * length(input$"_readout_gene")
  })
  
  # Data
  num_reads_thr <- 50  # threshold for UMI count in no-UMI assays
  min_dist_to_median <- 0.0
  
  observeEvent(input$file, {
    
    infile <- input$file
    if (is.null(infile)){
      v$data <- NULL    
    }
    df <- read_csv(infile$datapath, col_types = "cccccccdici")
    df <- df %>% 
      mutate(replicate = paste0("rep", replicate)) %>% 
      mutate(day = paste("day", day, sep = "_")) %>% 
      mutate(measurement = paste(replicate, assay, sep = "_")) %>% 
      mutate(plate = paste(measurement, day, sep = "_")) %>% 
      mutate(measurement = factor(measurement))
    
    # Filtering for assay quality
    
    df <- df %>% filter((assay == "no-UMI" & num_reads >= num_reads_thr) | assay == "UMI")
    df <- df %>% filter(!grepl("LOW|VL|NC", cell_density))
    
    v$data <- df
  })
  
  # Plot
  
  plot <- eventReactive(input$submit, {
    # Gather data & parameters
    df <- v$data
    pert_ <- input$"_pert"
    pert_type_ <- input$"_pert_type"
    rep_ <- input$"_rep"
    day_ <- input$"_day"
    readout_gene_ <- input$"_readout_gene"
    
    if (is.null(df)) return(NULL)
    if (is.null(pert_)) return(NULL)
    if (is.null(pert_type_)) return(NULL)
    if (is.null(rep_)) return(NULL)
    if (is.null(day_)) return(NULL)
    if (is.null(readout_gene_)) return(NULL)
    
    # Filter master df by selected parameters
    df <- df %>% 
      filter(  target_gene == pert_ &
                 shrna %in% pert_type_ &
                 replicate %in% gsub("(.+?)(\\_.*)", "\\1", rep_) &
                 assay %in% gsub("^.*\\_","",rep_) &
                 day %in% day_ &
                 readout_gene %in% readout_gene_)
    
    make_summary_plot(df, pert_, "MAE")
  })
  
  output$mainplotui <- renderUI({
      plotOutput("mainplot",
                 width = tab_width(),
                 height = tab_height())
  })
  
  output$mainplot <- renderPlot({
    plot()
  })
  
  # Reactively get data to populate UI
  
  ## Independent
  getPert <- reactive({
    df <- v$data
    if (is.null(df)) return(NULL)
    unique(df$target_gene)
  })
  getDay <- reactive({
    df <- v$data
    if (is.null(df)) return(NULL)
    # "day_{n}" --> n
    # as.numeric(gsub("^.*\\_","",unique(df$day)))
    unique(df$day)
  })
  
  ## Dependent on Perturbation selected
  getPertType <- reactive({
    df <- v$data
    pert <- input$"_pert"
    if (is.null(df)) return(NULL)
    if (is.null(pert)) return(NULL)
    
    shrnas <- df %>% 
      filter(target_gene == pert) %>%
      pull(shrna)
    
    unique(shrnas)
  })
  getRep <- reactive({
    df <- v$data
    pert <- input$"_pert"
    if (is.null(df)) return(NULL)
    if (is.null(pert)) return(NULL)
    
    rep_assays <- df %>% 
      filter(target_gene == pert) %>%
      mutate(rep_assay=paste0(replicate,"_",assay)) %>% 
      pull(rep_assay)
    
    unique(rep_assays)
  })
  getReadOutGene <- reactive({
    df <- v$data
    pert <- input$"_pert"
    if (is.null(df)) return(NULL)
    if (is.null(pert)) return(NULL)
    
    readout_genes <- df %>% 
      filter(target_gene == pert) %>%
      pull(readout_gene)
    
    unique(readout_genes)
  })
  
  # Render UI for parameters
  
  output$pert <- renderUI({
    selectInput('_pert', 'Perturbation', getPert())
  })
  output$pert_type <- renderUI({
    pickerInput('_pert_type', 'Perturbation type', 
                choices = getPertType(),
                selected = getPertType(),
                options = list(
                  `actions-box` = TRUE, 
                  size = 10,
                  `selected-text-format` = "count > 3"
                ), 
                multiple = TRUE)
  })
  output$rep <- renderUI({
    pickerInput('_rep', 'Replicant', 
                choices = getRep(),
                selected = getRep(),
                options = list(
                  `actions-box` = TRUE, 
                  size = 10,
                  `selected-text-format` = "count > 3"
                ), 
                multiple = TRUE)
  })
  output$day <- renderUI({
    pickerInput('_day', 'Day', 
                choices = getDay(),
                selected = getDay(),
                options = list(
                  `actions-box` = TRUE, 
                  size = 10
                ), 
                multiple = TRUE)
  })
  output$readout_gene <- renderUI({
    pickerInput('_readout_gene', 'Readout Gene', 
                choices = getReadOutGene(),
                selected = getReadOutGene(),
                options = list(
                  `actions-box` = TRUE, 
                  size = 10,
                  `selected-text-format` = "count > 3"
                ), 
                multiple = TRUE)
  })
  
  # Download plot
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$dataset, '.png', sep='') },
    content = function(file) {
      ggsave(file, plot = plot(), device = "png", 
             width = tab_width()/100, height = tab_height()/100, units = "in", limitsize=FALSE)
    }
  )
})

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
  
  # Plot dimensions
  ## Main screen plot
  width_per_day = 250
  height_per_replicate = 20
  ## By drug plot 
  width_of_readout_genes = 500
  height_per_day = 250
  
  # Store reactive values
  v <- reactiveValues(
    data = NULL,
    makeplot = FALSE
  )
  
  tab_width = eventReactive(input$submit, {
    
    width_per_day * nrow(v$data %>% filter(day %in% input$"_day") %>% count(day))
    # width_per_day * length(unique(v$filtered_df$day))
    # width_per_day * length(input$"_day")
  })
  tab_height = eventReactive(input$submit, {
    
    len_of_rep <- length(input$"_rep")
    len_of_readout <- nrow(v$data %>% filter(readout_gene %in% input$"_readout_gene") %>% count(readout_gene))
    
    len_of_rep * height_per_replicate * len_of_readout
    
    # row_height <- length(unique(v$filtered_df$replicate)) * height_per_replicant
    # row_height * length(unique(v$filtered_df$readout_gene))
    # 
    # row_height <- length(input$"_rep") * height_per_replicant
    # row_height * length(input$"_readout_gene")
  })
  drug_tab_height = eventReactive(input$submit_drug, {
    
    height_per_day * nrow(v$data %>% filter(day %in% input$"_day_drug") %>% count(day))

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
  
  ## Main Plot
  
  check_data <- function(data) {
    
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
    
    return (df)
  }
  
  plot <- eventReactive(input$submit, {
    
    df <- check_data(v$data)
    df <- df %>% 
      filter(  target_gene == input$"_pert" &
                 shrna %in% input$"_pert_type" &
                 replicate %in% gsub("(.+?)(\\_.*)", "\\1", input$"_rep") &
                 assay %in% gsub("^.*\\_","",input$"_rep") &
                 day %in% input$"_day" &
                 readout_gene %in% input$"_readout_gene")
    if (is.null(df)) return (NULL)
    
    switch (input$plottype,
      "ridge" = make_ridge_plot(df, input$"_pert", "MAE"),
      "violin" = make_violin_plot(df, input$"_pert", "MAE"),
      "box" = make_box_plot(df, input$"_pert", "MAE")
    )
    
  })
  
  output$mainplot <- renderPlot({
    plot()
  })
  
  output$mainplotui <- renderUI({
    
    if (tab_width() == 0 || tab_height() == 0) {
      return(NULL)
    }
    
      plotOutput("mainplot",
                 width = tab_width(),
                 height = tab_height())
  })
  
  ## Drug Plot
  
  drug_plot <- eventReactive(input$submit_drug, {
    
    df <- check_data(v$data)
    df <- df %>% 
      filter(day %in% input$"_day_drug" & readout_gene %in% input$"_readout_gene_drug")
    if (is.null(df)) return (NULL)
    
    make_drug_plot(df, input$"_readout_gene_drug", "MAE")
    
  })
  
  output$drugplot <- renderPlot({
   drug_plot()
  })
  
  output$drugplotui <- renderUI({
    
    if (tab_width() == 0 || tab_height() == 0) {
      return(NULL)
    }
    
    plotOutput("drugplot",
               width = width_of_readout_genes,
               height = drug_tab_height())
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
  
  ## Main plot parameters
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
    pickerInput('_rep', 'Replicate', 
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
  
  ## Drug plot parameters
  
  output$readout_gene_drug <- renderUI({
    selectInput('_readout_gene_drug', 'Perturbation', getReadOutGene())
  })
  output$day_drug <- renderUI({
    pickerInput('_day_drug', 'Day', 
                choices = getDay(),
                selected = getDay(),
                options = list(
                  `actions-box` = TRUE, 
                  size = 10
                ), 
                multiple = TRUE)
  })
  
  # Download plot
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0(tools::file_path_sans_ext(input$file$datapath), '.png') },
    content = function(file) {
      ggsave(file, plot = plot(), device = "png", 
             width = tab_width()/100, height = tab_height()/100, units = "in", limitsize=FALSE)
    }
  )
  
  output$downloadDrugPlot <- downloadHandler(
    filename = function() { paste0(tools::file_path_sans_ext(input$file$datapath), '.drug.png') },
    content = function(file) {
      ggsave(file, plot = drug_plot(), device = "png", 
             width = width_of_readout_genes/100, height = drug_tab_height()/100, units = "in", limitsize=FALSE)
    }
  )
})

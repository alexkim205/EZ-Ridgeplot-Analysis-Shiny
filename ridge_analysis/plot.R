# Taken from make_screen_summary_ridgeplot.R
# Slightly modified for R Shiny use

#! Plots bias distribution for all perturbations except select_target_gene overlayed with bias for shRNA against select_target_gene, 
#! divided by readout gene, measurement, and day of measurement.
#! @param df Dataframe with columns readout_class, readout_gene, measurement (should be of class "factor"), day, target_gene, shrna, bias
#! @param select_readout_class Readout class to display, e.g. "MAE", "BAE"
#! @output Plot for select_target_gene

plot_theme <- theme(text = element_text(size = 15), 
                    strip.text.y = element_text(size = 12, angle = 270), 
                    strip.text = element_text(face="bold", size=15, margin = margin(.2, .2, .2, .2, "cm")),
                    panel.grid.major = element_line(colour = "white",size=0.2), 
                    panel.spacing.x = unit(1, "lines"),
                    # axis.text.y = element_text(size = 12, angle = 30, vjust = -1.8),
                    axis.text.x = element_text(size = 12))

make_violin_plot <- function(df, select_target_gene, select_shrna, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% 
    filter(target_gene == select_target_gene, shrna %in% select_shrna)
  
  # Plotting
  
  plot <- ggplot(df, aes(measurement, bias)) + 
    facet_grid(readout_gene ~ day) + 
    labs(x = "measurement", y = "bias", title = select_target_gene) +
    
    geom_violin() + coord_flip() + 
    geom_segment(data = df_target, size = 1, aes(y = bias, yend = bias, x = as.numeric(measurement)-0.3, xend = as.numeric(measurement) + 0.3, color = shrna)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks=c(0.0,0.25,0.50,0.75,1.00), labels=c("0","0.25","0.50","0.75","1")) +
    plot_theme
  
  return(plot)
}

make_ridge_plot <- function(df, select_target_gene, select_shrna, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% 
    filter(target_gene == select_target_gene, shrna %in% select_shrna)
  
  # Plotting
  
  plot <- ggplot(df, aes(bias, measurement)) + 
    facet_grid(readout_gene ~ day) + 
    labs(x = "bias", y = "measurement", title = select_target_gene) +
    geom_density_ridges(scale = 0.9, size = 0.2) +
    geom_segment(data = df_target, aes(x = bias, xend = bias, y = as.numeric(measurement), yend = as.numeric(measurement) + 0.9, color = shrna)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks=c(0.0,0.25,0.50,0.75,1.00), labels=c("0","0.25","0.50","0.75","1")) +
    plot_theme
  
  return(plot)
}

make_box_plot <- function(df, select_target_gene, select_shrna, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% 
    filter(target_gene == select_target_gene, shrna %in% select_shrna)
  
  # Plotting
  
  plot <- ggplot(df, aes(measurement, bias)) + 
    facet_grid(readout_gene ~ day) + 
    labs(x = "measurement", y = "bias", title = select_target_gene) +
    geom_boxplot(alpha = 0.5, color = '#696969') + coord_flip() + 
    geom_segment(data = df_target, size = 1, aes(y = bias, yend = bias, x = as.numeric(measurement)-0.3, xend = as.numeric(measurement) + 0.3, color = shrna)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks=c(0.0,0.25,0.50,0.75,1.00), labels=c("0","0.25","0.50","0.75","1")) +
    plot_theme
  
  return(plot)
}

make_drug_plot <- function(df, select_target_gene, select_shrna, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% 
    filter(target_gene == select_target_gene, shrna %in% select_shrna)
  
  ggplot(df, aes(target_gene, bias)) + 
    facet_grid(day ~ readout_gene) + 
    geom_point() + 
    geom_point(data = df_target, aes(color = shrna)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    plot_theme
}
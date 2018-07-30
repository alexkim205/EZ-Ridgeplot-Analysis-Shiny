library(tidyverse)
library(ggridges)

# Variables

num_reads_thr <- 50  # threshold for UMI count in no-UMI assays
min_dist_to_median <- 0.0

# Data munging
path <- "/Users/alexkim/Dropbox/Gimelbrant_Lab/ridgeplots/by_alex"
file2 <- "data/20180403_Drug_screen_summary.csv"
file <- "data/20180306_DBI31_sequencing_summary.csv"
#df <- read_csv(paste0(path, "data/test_data.csv"), col_types = "cccccccdii")
df <- read_csv(file.path(path, "data/20180306_DBI31_sequencing_summary.csv"), col_types = "cccccccdici")

df <- df %>% 
  mutate(replicate = paste0("rep", replicate)) %>% 
  mutate(day = paste("day", day, sep = "_")) %>% 
  mutate(measurement = paste(replicate, assay, sep = "_")) %>% 
  mutate(plate = paste(measurement, day, sep = "_")) %>% 
  mutate(measurement = factor(measurement))

# Filtering for assay quality

df <- df %>% filter((assay == "no-UMI" & num_reads >= num_reads_thr) | assay == "UMI")
df <- df %>% filter(!grepl("LOW|VL|NC", cell_density))

plot_theme <- theme(text = element_text(size = 15), 
                    strip.text.y = element_text(size = 12, angle = 270), 
                    strip.text = element_text(face="bold", size=15, margin = margin(.2, .2, .2, .2, "cm")),
                    panel.grid.major = element_line(colour = "white",size=0.2), 
                    # axis.text.y = element_text(size = 12, angle = 30, vjust = -1.8),
                    axis.text.x = element_text(size = 12))

# Plot

make_violin_plot <- function(df, select_target_gene, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% filter(target_gene == select_target_gene)
    
  # Plotting
  
  plot <- ggplot(df, aes(measurement, bias)) + 
    facet_grid(readout_gene ~ day) + 
    labs(x = "measurement", y = "bias", title = select_target_gene) +
    # geom_density_ridges(scale = 0.9, size = 0.2) +
    
    geom_violin() + coord_flip() + 
    # geom_jitter(alpha=0.1) + 
    # geom_boxplot(aes(x=measurement, fill=measurement), alpha=0, lwd=0.2) + coord_flip() +
    geom_segment(data = df_target, size = 1, aes(y = bias, yend = bias, x = as.numeric(measurement)-0.3, xend = as.numeric(measurement) + 0.3, color = shrna)) + 
    # geom_point(data=df_target, aes(y=bias, x=measurement, color=shrna), size=2, alpha=0.7) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks=c(0.0,0.25,0.50,0.75,1.00), labels=c("0","0.25","0.50","0.75","1")) +
    plot_theme
  
  return(plot)
}

make_ridge_plot <- function(df, select_target_gene, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% filter(target_gene == select_target_gene)
  
  # Plotting
  
  plot <- ggplot(df, aes(bias, measurement)) + 
    facet_grid(readout_gene ~ day) + 
    labs(x = "bias", y = "measurement", title = select_target_gene) +
    geom_density_ridges(scale = 0.9, size = 0.2) +
    geom_segment(data = df_target, size = 1, aes(x = bias, xend = bias, y = as.numeric(measurement), yend = as.numeric(measurement) + 0.9, color = shrna)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks=c(0.0,0.25,0.50,0.75,1.00), labels=c("0","0.25","0.50","0.75","1")) +
    plot_theme
  
  return(plot)
}

make_box_plot <- function(df, select_target_gene, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% filter(target_gene == select_target_gene)
  
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

make_drug_plot <- function(df, select_target_gene, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day))))) %>%
    mutate(bias = as.numeric(bias))
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(target_gene!= select_target_gene)
  df_target <- df %>% filter(target_gene == select_target_gene)
  
  ggplot(df, aes(target_gene, bias)) + 
    facet_grid(day ~ readout_gene) + 
    geom_point() + 
    geom_point(data = df_target, aes(color = shrna)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    plot_theme
}

make_drug_plot(df %>% filter(readout_gene == "Dnajc12"), "Smchd1", "MAE")

make_box_plot(df %>% filter(readout_gene == "Adamtsl4"), "Smchd1", "MAE")
make_violin_plot(df %>% filter(readout_gene == "Adamtsl4"), "Smchd1", "MAE")
make_ridge_plot(df %>% filter(readout_gene == "Adamtsl4"), "Smchd1", "MAE")

make_summary_plot(df, "Dnmt1", "MAE")

target_genes <- unique(df$target_gene)
for (select_target_gene in target_genes) {
  make_summary_plot(df, select_target_gene, "MAE")
  }

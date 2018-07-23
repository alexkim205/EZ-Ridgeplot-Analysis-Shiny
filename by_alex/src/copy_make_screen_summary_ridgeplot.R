library(tidyverse)
library(ggridges)

# Variables

num_reads_thr <- 50  # threshold for UMI count in no-UMI assays
min_dist_to_median <- 0.0

# Data munging

# df <- read_csv("../data/test_data.csv", col_types = "cccccccdii")
df <- read_csv("../data/20180404_Drug_screen_summary.csv", col_types = "cccccccdici") # 
df <- df %>% 
  mutate(replicate = paste0("rep", replicate)) %>% 
  mutate(day = paste("day", day, sep = "_")) %>% 
  mutate(measurement = paste(replicate, assay, sep = "_")) %>% 
  mutate(plate = paste(measurement, day, sep = "_")) %>% 
  mutate(measurement = factor(measurement))

# Filtering for assay quality

df <- df %>% filter((assay == "no-UMI" & num_reads >= num_reads_thr) | assay == "UMI")
df <- df %>% filter(!grepl("DIFFUSE", cell_density))

# Plot

make_summary_plot <- function(df, select_Drug, select_readout_class = NA) {
  #! Plots bias distribution for all perturbations except select_Drug overlayed with bias for Concentration against select_Drug, 
  #! divided by readout gene, measurement, and day of measurement.
  #! @param df Dataframe with columns readout_class, readout_gene, measurement (should be of class "factor"), day, Drug, Concentration, bias
  #! @param select_readout_class Readout class to display, e.g. "MAE", "BAE"
  #! @output Plot for select_Drug
  
  # Filtering for genes of interest
  if (!is.na(select_readout_class)) {
    df <- df %>% filter(readout_class == select_readout_class)
  }
  df_non_target <- df %>% filter(Drug!= select_Drug)
  df_target <- df %>% filter(Drug == select_Drug)

  # Plotting

  plot <- ggplot(df_non_target, aes(bias, measurement)) + 
    facet_grid(readout_gene ~ day) + 
    labs(x = "bias",
         y = "measurement",
         title = select_Drug
    ) +
    # geom_density_ridges(aes(height=..density..), scale = 0.9, stat = "density", size = 0.2) +
    geom_density_ridges(scale = 0.9, size = 0.2) +
    geom_segment(data = df_target, aes(x = bias, xend = bias, y = as.numeric(measurement), yend = as.numeric(measurement) + 0.9, color = Concentration)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
    theme(text = element_text(size = 5), strip.text.y = element_text(angle = 0), 
          panel.grid.major =   element_line(colour = "white",size=0.2), 
          panel.spacing.x = unit(1, "lines"))
  
  pdf(paste0("../plots/", select_Drug, "_ridgeplot.pdf"))
  print(plot)
  dev.off()
  
}

# make_summary_plot(df, "Smchd1", "MAE")
# make_summary_plot(df, "Dnmt1", "MAE")

Drugs <- unique(df$Drug)
for (select_Drug in Drugs) {make_summary_plot(df, select_Drug, "MAE")}

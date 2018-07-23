library(tidyverse)
library(ggridges)

# Variables

num_reads_thr <- 50  # threshold for UMI count in no-UMI assays
min_dist_to_median <- 0.0

# Data munging
path <- "/Users/alexkim/Desktop/Gimelbrant/ridgeplots/by_alex/"
file <- "data/20180306_DBI31_sequencing_summary.csv"
#df <- read_csv(paste0(path, "data/test_data.csv"), col_types = "cccccccdii")
df <- read_csv(paste0(path, "data/20180306_DBI31_sequencing_summary.csv"), col_types = "cccccccdici")

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

make_summary_plot <- function(df, select_target_gene, select_readout_class = NA) {
  
  # Prettify
  df <- df %>% mutate(day = paste0("Day ", sprintf("%02d",as.numeric(gsub("^.*\\_","",day)))))
  
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
    geom_jitter(alpha=0.1) + 
    # geom_boxplot(aes(x=measurement, fill=measurement), alpha=0, lwd=0.2) + coord_flip() +
    geom_segment(data = df_target, aes(y = bias, yend = bias, x = as.numeric(measurement)-0.3, xend = as.numeric(measurement) + 0.3, color = shrna)) + 
    # geom_point(data=df_target, aes(y=bias, x=measurement, color=shrna), size=2, alpha=0.7) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks=c(0.0,0.25,0.50,0.75,1.00), labels=c("0","0.25","0.50","0.75","1")) +
    plot_theme
  
  return(plot)
}

make_summary_plot(df %>% filter(readout_gene == "Adamtsl4"), "Smchd1", "MAE")
make_summary_plot(df, "Dnmt1", "MAE")

target_genes <- unique(df$target_gene)
for (select_target_gene in target_genes) {
  make_summary_plot(df, select_target_gene, "MAE")
  }

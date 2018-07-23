library(tidyverse)

# Variables

num_reads_thr <- 50  # threshold for UMI count in no-UMI assays
min_dist_to_median <- 0.0

# Data munging

# df <- read_csv("../data/test_data.csv", col_types = "cccccccdii")
df <- read_csv("../data/20180306_DBI31_sequencing_summary.csv", col_types = "cccccccdici")
df <- df %>% mutate(replicate = paste0("rep", replicate))
df <- df %>% mutate(day = paste("day", day, sep = "_"))
df <- df %>% mutate(measurement = paste(replicate, assay, sep = "_"))
df <- df %>% mutate(plate = paste(measurement, day, sep = "_"))

# Filtering for assay quality

df <- df %>% filter((assay == "no-UMI" & num_reads >= num_reads_thr) | assay == "UMI")
df <- df %>% filter(!grepl("LOW|VL|NC", cell_density))

# Double MAD normalization of bias
# For explanation, see https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

df <- df %>% group_by(plate, readout_gene) %>% 
  mutate(abs_bias = abs(bias - 0.5)) %>% 
  mutate(abs_bias_median = median(abs_bias, na.rm = TRUE)) %>% 
  mutate(abs_bias_mad = median(abs(abs_bias - abs_bias_median)), na.rm = TRUE) %>% 
  mutate(abs_bias_halfmad = median(abs(abs_bias[abs_bias < abs_bias_median] - abs_bias_median[abs_bias < abs_bias_median]), na.rm = TRUE)) %>% 
  ungroup
df <- df %>% mutate(halfmad_score = ifelse(abs_bias < abs_bias_median, abs(abs_bias - abs_bias_median) / abs_bias_halfmad, 0))
df <- df %>% mutate(lom_confidence = ifelse(abs_bias < abs_bias_median - min_dist_to_median, abs(abs_bias - abs_bias_median) / abs_bias_halfmad, 0))  # requires difference to median to be at least min_dist_to_median

# Plot

make_summary_plot <- function(df, select_target_gene, select_readout_class) {
  
  # Filtering for genes of interest
  
  df_target_gene <- df %>% filter(readout_class == select_readout_class)
  df_target_gene <- df %>% filter(target_gene == select_target_gene)
  
  # Plotting

  plot <- ggplot(df_target_gene, aes(shrna, measurement, fill = lom_confidence)) + 
    geom_tile(colour = "white") + 
    facet_grid(readout_gene ~ day) + 
    scale_fill_gradient(low = "white", high = "green", trans = "log1p", limits = c(0, 20), na.value = "green") +
    labs(x = "shRNA",
         y = "measurement",
         title = select_target_gene
    ) + 
    theme(text = element_text(size = 5), strip.text.y = element_text(angle = 0))
  
  pdf(paste0("../plots/", select_target_gene, ".pdf"))
  print(plot)
  dev.off()
  
}

# make_summary_plot(df, "Smchd1", "MAE")
# make_summary_plot(df, "Dnmt1", "MAE")

target_genes <- unique(df$target_gene)
for (select_target_gene in target_genes) {make_summary_plot(df, select_target_gene, "MAE")}

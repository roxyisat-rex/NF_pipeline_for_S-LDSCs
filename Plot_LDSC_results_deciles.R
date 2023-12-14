# Load required libraries
library(ggplot2)
library(cowplot)
library(dplyr)

# Function to plot LDSC results
plot_ldsc_results <- function(input_csv, output_pdf) {
  # Read the LDSC results CSV file
  ldsc_results <- read.csv(input_csv)

  # Assuming you have a column 'Minus_Log10_P' containing the calculated -log10(p) values
  # Adjust this based on the actual column name in your dataset
  ldsc_results$Quantile <- as.character(ldsc_results$Quantile)  # Ensure Quantile is treated as a character

  # Plot the results
  p <- ggplot(ldsc_results, aes(x = Quantile, y = Minus_Log10_P)) +
    geom_point() +
    geom_line() +
    facet_wrap(~Cell_Type, scales = "free_y", labeller = labeller(Cell_Type = label_both)) +  # Display cell type names
    theme_cowplot(12, strip.background = element_rect(fill = "#89CFF0")) +  # Use theme_cowplot with font size 12 and strip background color
    labs(x = "Quantile", y = "-Log10(P)", title = "LDSC Results for Each Cell Type") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  # Save the plot as a PDF
  ggsave(output_pdf, plot = p, width = 12, height = 8)
}

# Run the function with command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  plot_ldsc_results(args[1], args[2])
} else {
  cat("Usage: Rscript plot_ldsc_results.R <input_csv_file> <output_pdf_file>\n")
}

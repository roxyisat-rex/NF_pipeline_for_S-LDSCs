#!/usr/bin/env Rscript

library(FSA)

# Check if command-line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) < 2) {
  stop("Usage: Rscript NF_multiple_corrections.R <input_file> <output_file>")
}

# Get input and output file paths from command-line arguments
input_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

Multi_corrections <- function(input_file, output_file) {
  data <- read.csv(input_file)
  df_data <- data.frame(data)
  corrected_data <- p.adjust(df_data$p_value, method = "bonferroni")
  results <- cbind(df_data, corrected_data)
  colnames(results)[6] <- "P.Adj"
  minus.log10 <- -log10(results$P.Adj)
  results$minus_log10 <- minus.log10
  colnames(results)[1] <- "Cell_types"
  write.csv(results, row.names = FALSE, output_file)
}

# Call the function with provided input and output file paths
Multi_corrections(input_file, output_file)

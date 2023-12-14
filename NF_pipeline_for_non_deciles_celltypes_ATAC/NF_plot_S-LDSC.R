# plot_analysis.R

plot_analysis <- function(input_csv, output_pdf) {
  library(ggplot2)
  library(ggpubr)
  library(grid)
  library(patchwork)
  library(cowplot)
  
  # Read data from CSV
  df <- read.csv(input_csv)
  
  #PLOT the -log 10 
  df.p<- ggplot(df, aes(y=(minus_log10), x=Cell_types)) + geom_bar(aes(fill = Cell_types), stat = 'identity') + scale_x_discrete(expand = c(0,0)) + ylab("-Log 10") + xlab(" ") + scale_y_continuous(expand = c(0,0)) + labs(title = " ")
  final.p <- df.p + theme_cowplot(12) + theme(axis.text.x = element_text(angle = 70, vjust =1, hjust = 1), strip.background = element_rect(fill="#89CFF0"), plot.title = element_text(hjust = 0.5), legend.position = "none") + guides(fill = guide_legend(title = "Cell types ")) + geom_hline(yintercept = 1.301029996, linetype = "dashed", color = "red")
  final.p
  
  #plot the enrichment 
  df.p1<- ggplot(df, aes(y=(Enrichment), x=Cell_types)) + geom_bar(aes(fill = Cell_types), stat = 'identity') + scale_x_discrete(expand = c(0,0)) + ylab("Enrichment") + xlab(" ") + scale_y_continuous(expand = c(0,0)) + labs(title = " ") 
  final.p1 <- df.p1 + theme_cowplot(12) + theme(axis.text.x = element_text(angle = 70, vjust =1, hjust = 1), strip.background = element_rect(fill="#89CFF0"), plot.title = element_text(hjust = 0.5), legend.position = "none") + guides(fill = guide_legend(title = "Cell types ")) 
  final.p1
  
  figure <- ggarrange (final.p,final.p1, ncol = 2) #, common.legend = T, legend = "bottom")
  figure <- annotate_figure(figure, top = text_grob("S-LDSC results", face = "bold", size = 16))
  # Save the final figure
  graphics.off()
  ggsave(output_pdf, figure, width = 12, height = 8)
}

# Run the function with command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  plot_analysis(args[1], args[2])
} else {
  cat("Usage: Rscript plot_analysis.R <input_csv_file> <output_pdf_file>\n")
}

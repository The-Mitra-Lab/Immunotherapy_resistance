# load necessary libraries
library(ggplot2)

# get the input file from the command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input file provided. Usage: Rscript plot_insert_size.R <input_file>")
}
input_file <- args[1]

# read data
data <- read.table(input_file, header = FALSE, col.names = c("Count", "FragmentLength"))

# filter data for fragment lengths between 0 and 1000 bp
data <- subset(data, FragmentLength >= 0 & FragmentLength <= 1000)

# plot
p <- ggplot(data, aes(x = FragmentLength, y = Count)) +
  geom_line(color = "blue") +  # Line plot
  geom_vline(xintercept = 150, linetype = "dashed", color = "red", size = 1) +  # Vertical line at x = 150
  theme_minimal() +
  labs(title = "Fragment Length Distribution (0-1000 bp)",
       x = "Fragment Length (bp)",
       y = "Count")

# generate output file name
output_file <- paste0("./QC/", tools::file_path_sans_ext(basename(input_file)), "_fragment_length_plot.png")

# save
ggsave(output_file, plot = p, width = 6, height = 4, dpi = 300)

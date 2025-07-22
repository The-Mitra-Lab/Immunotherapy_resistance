args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_deseq2.R <combinedBed>")
}

# define input and output file names
input_file <- args[1]
output_all <- "deseq2_results_all.tsv"
output_sig_up <- "deseq2_results_up_significant.tsv"
output_sig_down <- "deseq2_results_down_significant.tsv"

# load required libraries
library(DESeq2)
library(tidyverse)

# read in merged peak count matrix with custom column names
cols <- c("chr", "start", "end", "wt1", "wt2", "ko1", "ko2")
df <- read_tsv(input_file, col_names = cols)

# create unique peak IDs and set as rownames
df <- df %>%
  mutate(peak_id = paste(chr, start, end, sep = "_")) %>%
  column_to_rownames("peak_id")

# extract just the count data
count_data <- df %>% select(wt1, wt2, ko1, ko2)

# create metadata for DESeq2 with condition labels
col_data <- data.frame(
  row.names = colnames(count_data),
  condition = c("wt", "wt", "ko", "ko")
)
col_data$condition <- factor(col_data$condition, levels = c("wt", "ko"))

# build and run DESeq2 model
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(count_data)),
                              colData = col_data,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# save all results
res_all <- as.data.frame(res) %>% rownames_to_column("peak_id")
write_tsv(res_all, output_all)

# save significantly upregulated peaks (padj < 0.05 and log2FC > 1)
res_sig <- res_all %>%
  filter(padj < 0.05, log2FoldChange > 1)
write_tsv(res_sig, output_sig_up)

# save significantly downregulated peaks (padj < 0.05 and log2FC < -1)
res_sig <- res_all %>%
  filter(padj < 0.05, log2FoldChange < -1)
write_tsv(res_sig, output_sig_down)

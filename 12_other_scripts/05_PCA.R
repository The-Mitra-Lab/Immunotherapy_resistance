# load libraries
library(DESeq2)
library(ggplot2)

# load count matrix 
counts <- read.table("path/to/union_peak_counts.txt", header = FALSE)

# assign sample names
sample_names <- c(
  "MOC1veh_restim_1", "MOC1veh_restim_2", "MOC1veh_1", "MOC1veh_2",
  "MOC1ifn_restim_1", "MOC1ifn_restim_2", "MOC1ifn_1", "MOC1ifn_2",
  "KOveh_restim_2", "KOveh_1", "KOveh_2", "KOifn_restim_2",
  "DS_KOifn_restim_1", "DS_KOveh_restim_1", "KOifn_1", "KOifn_2"
)
colnames(counts) <- c("chr", "start", "end", sample_names)

# set rownames as unique peak IDs
rownames(counts) <- make.unique(paste0(counts$chr, ":", counts$start, "-", counts$end))

# drop coordinate columns
counts <- counts[, sample_names]

# remove restim samples 
non_restim_samples <- sample_names[!grepl("restim", sample_names)]
counts_filtered <- counts[, non_restim_samples]

# define conditions
conditions <- factor(c(
  "MOC1veh", "MOC1veh",
  "MOC1ifn", "MOC1ifn",
  "KOveh", "KOveh",
  "KOifn", "KOifn"
))
col_data <- data.frame(row.names = non_restim_samples, condition = conditions)

# ----- DESeq2 workflow -----
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = col_data, design = ~condition)

# sum reads across all samples for each peak
peak_sums <- rowSums(counts(dds))

# histogram
hist(peak_sums, breaks = 100, col = "gray70",
     xlab = "Total reads per peak", main = "Distribution of peak read counts")
abline(v = 100, col = "red", lwd = 2, lty = 2)

# filter low-count peaks
keep <- rowSums(counts(dds) >= 100) >= 2
dds <- dds[keep, ]

# normalize with VST
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# keep top 1000 most variable peaks
vars <- apply(vst_mat, 1, var)
top_peaks <- order(vars, decreasing = TRUE)[1:min(1000, nrow(vst_mat))]
vst_top <- vst_mat[top_peaks, ]

# PCA 
pca <- prcomp(t(vst_top), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df$Condition <- col_data$condition[match(pca_df$Sample, rownames(col_data))]

# define custom colors
custom_colors <- c(
  "MOC1veh" = "black",
  "MOC1ifn" = "red",
  "KOveh"   = "royalblue",
  "KOifn"   = "green3"
)

# plot

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, shape = 21, stroke = 1.1, fill = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(values = custom_colors) +
  xlab(paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 2), "%)")) +
  ylab(paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 2), "%)")) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "none",
    plot.margin = margin(3, 3, 3, 3)
  )

# save plot
ggsave(
  filename = "path/to/pca_plot.png",
  width = 2.3,
  height = 2,
  dpi = 300
)

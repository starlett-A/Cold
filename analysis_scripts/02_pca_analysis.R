# =================================================================
# 📊 PCA ANALYSIS (UNSUPERVISED)
# Purpose: Pareto scaling + PCA on preprocessed metabolomics data
# Input:  Preprocessed data from 02_preprocessing_pipeline.R
# Output: PCA scores, explained variance, publication-quality figures
# =================================================================

cat("🚀 Starting PCA analysis...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

# Ensure pacman is available
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c(
  "dplyr", "tidyr", "readr", "ggplot2", "ggrepel",
  "mixOmics", "patchwork", "viridis", "magrittr"
))

cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD PREPROCESSED DATA
# =============================================
cat("\n📂 STEP 1: Loading preprocessed data...\n")

# Use consistent relative path (output from previous step)
input_file <- "output/preprocessed_data.csv"

if (!file.exists(input_file)) {
  stop("Error: Preprocessed data not found at: ", input_file,
       "\nPlease run 02_preprocessing_pipeline.R first.")
}

data <- read_csv(input_file, show_col_types = FALSE)

cat("Data loaded successfully:\n")
cat(sprintf("  Samples: %d\n", nrow(data)))
cat(sprintf("  Metabolites: %d\n", ncol(data) - 2))
cat(sprintf("  Groups: %s\n", paste(unique(data$Group), collapse = ", ")))

# =============================================
# 2. PARETO SCALING
# =============================================
cat("\n⚖️ STEP 2: Pareto scaling...\n")

# Extract numeric matrix
numeric_data <- as.matrix(data[, 3:ncol(data)])

# Pareto scaling function (x - mean) / sqrt(sd)
pareto_scale <- function(x) {
  x_mean <- mean(x, na.rm = TRUE)
  x_sd   <- sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(x - x_mean)   # fallback for zero-variance metabolites
  } else {
    return((x - x_mean) / sqrt(x_sd))
  }
}

scaled_numeric <- apply(numeric_data, 2, pareto_scale)

# Reassemble scaled data frame
scaled_data <- data.frame(
  Sample_ID = data$Sample_ID,
  Group     = data$Group,
  scaled_numeric,
  stringsAsFactors = FALSE
)

# Save scaled data for later use (e.g., PLS-DA)
write_csv(scaled_data, "output/preprocessed_data_scaled.csv")
cat("✅ Pareto scaling complete – saved to output/preprocessed_data_scaled.csv\n")

# =============================================
# 3. PCA COMPUTATION
# =============================================
cat("\n📊 STEP 3: PCA computation...\n")

pca_X <- as.matrix(scaled_data[, 3:ncol(scaled_data)])
rownames(pca_X) <- scaled_data$Sample_ID
pca_Y <- scaled_data$Group

pca_result <- mixOmics::pca(
  pca_X,
  ncomp  = 3,
  center = TRUE,   # already centered? harmless to center again
  scale  = FALSE   # data already Pareto scaled
)

# Variance explained
prop_var <- pca_result$prop_expl_var$X * 100
cat(sprintf("   PC1 explains %.1f%% variance\n", prop_var[1]))
cat(sprintf("   PC2 explains %.1f%% variance\n", prop_var[2]))

# Save proportion of variance
var_df <- data.frame(
  PC = paste0("PC", seq_along(prop_var)),
  Variance_Explained_Pct = prop_var
)
write_csv(var_df, "output/PCA_variance_explained.csv")
cat("✅ Variance explained saved: output/PCA_variance_explained.csv\n")

# =============================================
# 4. PCA SCORE PLOTS
# =============================================
cat("\n🎨 STEP 4: Generating PCA score plots...\n")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ---------- 4.1 Basic PCA plot (mixOmics style) ----------
pca_basic <- plotIndiv(
  pca_result,
  comp      = c(1, 2),
  group     = pca_Y,
  ind.names = FALSE,
  ellipse   = TRUE,
  legend    = TRUE,
  title     = "PCA Score Plot (Unsupervised)",
  X.label   = sprintf("PC1 (%.1f%%)", prop_var[1]),
  Y.label   = sprintf("PC2 (%.1f%%)", prop_var[2]),
  style     = "ggplot2"
)

ggsave("output/figures/PCA_Score_Plot_Basic.png",
       plot = pca_basic$graph, width = 10, height = 8, dpi = 600)
cat("✅ Basic PCA plot saved\n")

# ---------- 4.2 Publication-quality PCA plot ----------
cat("🎨 Creating publication-quality PCA plot...\n")

# Extract scores
pca_scores <- as.data.frame(pca_result$variates$X[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$Group   <- pca_Y
pca_scores$Sample  <- rownames(pca_scores)

# Group centroids (optional, kept for possible labeling)
group_centroids <- pca_scores %>%
  group_by(Group) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )

journal_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  # 95% confidence ellipses
  stat_ellipse(
    geom      = "polygon",
    alpha     = 0.1,
    level     = 0.95,
    linetype  = "dashed",
    linewidth = 0.3      # use linewidth instead of size for lines (ggplot2 >=3.4.0)
  ) +
  # Sample points
  geom_point(
    aes(shape = Group),
    size   = 4,
    alpha  = 0.8,
    stroke = 0.8
  ) +
  # Axis labels
  labs(
    x     = sprintf("PC1 (%.1f%%)", prop_var[1]),
    y     = sprintf("PC2 (%.1f%%)", prop_var[2]),
    title = "Principal Component Analysis Score Plot"
  ) +
  # Color and shape scales
  scale_color_brewer(palette = "Set2", name = "Group") +
  scale_fill_brewer(palette = "Set2", name = "Group") +
  scale_shape_manual(values = c(16, 17, 15, 18, 8)[1:length(unique(pca_Y))]) +
  # Theme (bold axis style)
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5,
                              margin = margin(b = 10)),
    axis.title   = element_text(size = 24, face = "bold", color = "black"),
    axis.text    = element_text(size = 18, face = "bold", color = "black"),
    axis.line    = element_line(color = "black", linewidth = 1.2),
    axis.ticks   = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position  = "top",
    legend.title     = element_text(size = 22),
    legend.text      = element_text(size = 20),
    legend.box       = "horizontal",
    legend.margin    = margin(t = 0, b = 10),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = unit(c(1, 1, 1, 1), "cm")
  ) +
  coord_fixed(ratio = 1)   # equal unit lengths on both axes

# Save as PNG and PDF
ggsave("output/figures/PCA_Score_Plot_Journal.png",
       plot = journal_pca, width = 10, height = 8, dpi = 600, bg = "white")
ggsave("output/figures/PCA_Score_Plot_Journal.pdf",
       plot = journal_pca, width = 10, height = 8, bg = "white")
cat("✅ Publication-quality PCA plot saved (PNG & PDF)\n")

# =============================================
# 5. FINAL SUMMARY
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("✅ PCA ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 50), "\n\n"))

cat("📊 OUTPUT FILES:\n")
cat("  - output/preprocessed_data_scaled.csv\n")
cat("  - output/PCA_variance_explained.csv\n")
cat("  - output/figures/PCA_Score_Plot_Basic.png\n")
cat("  - output/figures/PCA_Score_Plot_Journal.png\n")
cat("  - output/figures/PCA_Score_Plot_Journal.pdf\n")

cat("\n🔜 Next step: supervised analysis (PLS-DA)\n")
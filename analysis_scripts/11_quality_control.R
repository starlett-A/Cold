# =================================================================
# 🔬 QC SAMPLE ASSESSMENT (TECHNICAL REPRODUCIBILITY)
# Purpose: Evaluate QC sample clustering, detect outliers, and
#          compute RSD metrics to confirm data quality before
#          differential analysis.
# Input:   QC metabolite abundance file (data/QC_metabolites.csv)
# Output:  PCA plots, RSD distribution, quantitative assessment report
# =================================================================

cat("🚀 Starting QC sample assessment...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("tidyverse", "ggplot2", "ggrepel", "readr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD QC DATA
# =============================================
cat("\n📂 STEP 1: Loading QC sample data...\n")

qc_file <- "data/QC_metabolites.csv"
if (!file.exists(qc_file)) {
  stop("Error: QC data file not found at: ", qc_file,
       "\nPlease place your QC metabolite data in the 'data/' folder.")
}

qc_samples <- read_csv(qc_file, show_col_types = FALSE)

cat(sprintf("  QC samples : %d\n", nrow(qc_samples)))
cat(sprintf("  Columns    : %d\n", ncol(qc_samples)))

# Show a glimpse
cat("  First few columns:\n")
print(head(qc_samples[, 1:min(5, ncol(qc_samples))]))

# =============================================
# 2. PREPARE DATA FOR PCA
# =============================================
cat("\n🔧 STEP 2: Preparing numeric matrix for PCA...\n")

# Extract all numeric columns
pca_data <- qc_samples %>%
  select(where(is.numeric))

cat(sprintf("  Initial metabolite count: %d\n", ncol(pca_data)))

# Remove columns with any NA or zero variance
pca_data <- pca_data[, colSums(is.na(pca_data)) == 0]
pca_data <- pca_data[, apply(pca_data, 2, var) != 0]
cat(sprintf("  After cleaning: %d metabolites\n", ncol(pca_data)))

# Add sample IDs (if present)
if ("Sample_ID" %in% colnames(qc_samples)) {
  sample_ids <- qc_samples$Sample_ID
} else {
  sample_ids <- paste0("QC_", seq_len(nrow(pca_data)))
}

# =============================================
# 3. PERFORM PCA (scaled)
# =============================================
cat("\n📊 STEP 3: Running PCA (centered & scaled)...\n")

pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

# Extract scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample_ID <- sample_ids

# Variance explained
var_exp <- summary(pca_result)$importance[2, 1:2] * 100
cat(sprintf("  PC1 explains %.1f%% variance\n", var_exp[1]))
cat(sprintf("  PC2 explains %.1f%% variance\n", var_exp[2]))
cat(sprintf("  Cumulative (PC1+PC2): %.1f%%\n", var_exp[1] + var_exp[2]))

# =============================================
# 4. BASIC QC PCA PLOT
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 4: Generating QC PCA plot\n")
cat(paste0(strrep("=", 50), "\n"))

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, color = "#1E88E5", alpha = 0.7) +
  geom_text_repel(aes(label = Sample_ID), size = 3.2,
                  max.overlaps = Inf, box.padding = 0.5) +
  stat_ellipse(type = "norm", level = 0.95,
               color = "#D81B60", linetype = 2, linewidth = 0.8) +
  labs(
    title = "QC Sample Technical Assessment",
    subtitle = paste(nrow(pca_scores), "QC samples |", ncol(pca_data), "metabolites"),
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    caption = "Ideal: all QC samples tightly clustered, minimizing technical variation"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(pca_plot)

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
ggsave("output/figures/QC_PCA_basic.png", pca_plot, width = 10, height = 7, dpi = 600)
ggsave("output/figures/QC_PCA_basic.pdf", pca_plot, width = 10, height = 7)
cat("✅ Basic QC PCA plot saved\n")

# =============================================
# 5. QUANTITATIVE OUTLIER DETECTION
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 5: Outlier detection via distance to centre\n")
cat(paste0(strrep("=", 50), "\n"))

# Euclidean distance to origin (PC1=0, PC2=0)
pca_scores$dist_to_center <- sqrt(pca_scores$PC1^2 + pca_scores$PC2^2)
dist_threshold <- median(pca_scores$dist_to_center) + 2 * sd(pca_scores$dist_to_center)
pca_scores$is_outlier <- pca_scores$dist_to_center > dist_threshold

n_outliers <- sum(pca_scores$is_outlier)
cat(sprintf("  Distance threshold: %.2f\n", dist_threshold))
cat(sprintf("  Outliers detected : %d / %d\n", n_outliers, nrow(pca_scores)))

if (n_outliers > 0) {
  cat("  Outlier details:\n")
  print(pca_scores[pca_scores$is_outlier, c("Sample_ID", "PC1", "PC2", "dist_to_center")])
}

# =============================================
# 6. ENHANCED PCA PLOT WITH OUTLIERS
# =============================================
cat("\n🎨 Generating enhanced PCA plot (outliers highlighted)...\n")

pca_outlier_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = is_outlier, size = is_outlier), alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "#1E88E5", "TRUE" = "#D81B60"),
                     name = "Outlier?") +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5), guide = "none") +
  geom_text_repel(
    data = subset(pca_scores, is_outlier),
    aes(label = Sample_ID),
    size = 3.5, color = "#D81B60", box.padding = 0.8
  ) +
  stat_ellipse(type = "norm", level = 0.95, color = "darkgrey", linetype = 2) +
  labs(
    title = "PCA: QC Samples with Outliers Marked",
    subtitle = paste0("Outliers: ", n_outliers, " / ", nrow(pca_scores)),
    color = "Outlier?"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

print(pca_outlier_plot)
ggsave("output/figures/QC_PCA_outliers.png", pca_outlier_plot, width = 10, height = 7, dpi = 600)
ggsave("output/figures/QC_PCA_outliers.pdf", pca_outlier_plot, width = 10, height = 7)
cat("✅ Outlier PCA plot saved\n")

# =============================================
# 7. RSD CALCULATION (GOLD STANDARD)
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 6: RSD analysis (Relative Standard Deviation)\n")
cat(paste0(strrep("=", 50), "\n"))

# RSD = (sd / mean) * 100 for each metabolite
rsd_values <- apply(pca_data, 2, function(x) {
  if (length(x) < 3 || mean(x) == 0) return(NA)
  sd(x) / mean(x) * 100
})
rsd_values <- na.omit(rsd_values)

if (length(rsd_values) == 0) stop("No valid RSD values computed. Check your input data.")

cat(sprintf("  Metabolites analysed : %d\n", length(rsd_values)))
cat(sprintf("  Median RSD           : %.1f%%\n", median(rsd_values)))
cat(sprintf("  RSD < 20%%            : %.1f%%\n", sum(rsd_values < 20) / length(rsd_values) * 100))
cat(sprintf("  RSD < 30%%            : %.1f%%\n", sum(rsd_values < 30) / length(rsd_values) * 100))

# RSD distribution plot
rsd_plot <- ggplot(data.frame(RSD = rsd_values), aes(x = RSD)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of RSD in QC Samples",
    subtitle = sprintf("%.1f%% metabolites with RSD < 30%%",
                       sum(rsd_values < 30) / length(rsd_values) * 100),
    x = "Relative Standard Deviation (%)",
    y = "Number of Metabolites"
  ) +
  theme_minimal(base_size = 14)

ggsave("output/figures/QC_RSD_distribution.png", rsd_plot, width = 8, height = 6, dpi = 600)
ggsave("output/figures/QC_RSD_distribution.pdf", rsd_plot, width = 8, height = 6)
cat("✅ RSD distribution plot saved\n")

# =============================================
# 8. COMPREHENSIVE QUALITY REPORT
# =============================================
cat(paste0("\n", strrep("=", 70), "\n"))
cat("📋 COMPREHENSIVE QC QUALITY REPORT\n")
cat(paste0(strrep("=", 70), "\n"))

# PCA spread (SD instead of CV to avoid zero-mean issues)
cat("📊 PCA dispersion:\n")
cat(sprintf("  PC1 SD: %.2f (range %.2f to %.2f)\n",
            sd(pca_scores$PC1), min(pca_scores$PC1), max(pca_scores$PC1)))
cat(sprintf("  PC2 SD: %.2f (range %.2f to %.2f)\n",
            sd(pca_scores$PC2), min(pca_scores$PC2), max(pca_scores$PC2)))

# Distance statistics
cat("\n📏 Distance to centre:\n")
dist_stats <- summary(pca_scores$dist_to_center)
cat(sprintf("  Min    : %.1f\n", dist_stats["Min."]))
cat(sprintf("  Median : %.1f\n", dist_stats["Median"]))
cat(sprintf("  Mean   : %.1f\n", dist_stats["Mean"]))
cat(sprintf("  Max    : %.1f\n", dist_stats["Max."]))
cat(sprintf("  SD     : %.1f\n", sd(pca_scores$dist_to_center)))
dist_cv <- sd(pca_scores$dist_to_center) / mean(pca_scores$dist_to_center) * 100
cat(sprintf("  CV     : %.1f%%\n", dist_cv))

# Outlier summary
cat(sprintf("\n🔍 Outlier summary:\n"))
cat(sprintf("  Outliers         : %d/%d (%.1f%%)\n",
            n_outliers, nrow(pca_scores),
            n_outliers / nrow(pca_scores) * 100))
if (n_outliers > 0) {
  avg_dist <- mean(pca_scores$dist_to_center[!pca_scores$is_outlier])
  for (id in pca_scores$Sample_ID[pca_scores$is_outlier]) {
    dist_val <- pca_scores$dist_to_center[pca_scores$Sample_ID == id]
    cat(sprintf("  %s: distance = %.1f (%.1fx average of non-outliers)\n",
                id, dist_val, dist_val / avg_dist))
  }
}

# RSD evaluation
cat("\n🏆 RSD gold standard:\n")
cat(sprintf("  Metabolites RSD < 20%%: %.1f%%\n",
            sum(rsd_values < 20) / length(rsd_values) * 100))
cat(sprintf("  Metabolites RSD < 30%%: %.1f%%\n",
            sum(rsd_values < 30) / length(rsd_values) * 100))

# Overall assessment
cat("\n✅ Overall quality assessment:\n")
# RSD criteria (primary)
if (sum(rsd_values < 30) / length(rsd_values) >= 0.9) {
  cat("  RSD quality: Excellent (≥90% metabolites RSD<30%)\n")
} else if (sum(rsd_values < 30) / length(rsd_values) >= 0.7) {
  cat("  RSD quality: Good (≥70% metabolites RSD<30%)\n")
} else {
  cat("  RSD quality: Needs attention (<70% metabolites RSD<30%)\n")
}

# Outlier control
outlier_ratio <- n_outliers / nrow(pca_scores)
if (outlier_ratio <= 0.05) {
  cat("  Outlier control: Excellent (≤5% outliers)\n")
} else if (outlier_ratio <= 0.10) {
  cat("  Outlier control: Acceptable (5-10% outliers)\n")
} else {
  cat("  Outlier control: Needs attention (>10% outliers)\n")
}

# Distance CV
if (dist_cv < 50) {
  cat("  Spread uniformity: Good (distance CV < 50%)\n")
} else {
  cat("  Spread uniformity: High dispersion (distance CV ≥ 50%)\n")
}

# Final conclusion
cat("\n🎯 Final conclusion:\n")
if (sum(rsd_values < 30) / length(rsd_values) >= 0.9 && outlier_ratio <= 0.10) {
  cat("  Data quality is excellent. RSD meets strict criteria, only a few outliers.\n")
  cat("  Proceed with downstream differential analysis confidently.\n")
} else if (sum(rsd_values < 30) / length(rsd_values) >= 0.7) {
  cat("  Data quality is acceptable. RSD is satisfactory; consider checking outliers.\n")
} else {
  cat("  Data quality may need further preprocessing. Investigate QC outliers and RSD.\n")
}

cat(paste0(strrep("=", 70), "\n"))
cat("✅ QC ASSESSMENT COMPLETE\n")

# Save summary statistics as CSV
summary_list <- list(
  Metric = c("QC samples", "Metabolites after cleaning", "PC1 variance (%)", "PC2 variance (%)",
             "PC1 SD", "PC2 SD", "Median distance", "Mean distance", "Max distance", "Distance CV (%)",
             "Outliers (count)", "Outlier ratio", "Median RSD (%)", "RSD < 20% (%)", "RSD < 30% (%)"),
  Value = c(nrow(pca_scores), ncol(pca_data), var_exp[1], var_exp[2],
            sd(pca_scores$PC1), sd(pca_scores$PC2),
            median(pca_scores$dist_to_center), mean(pca_scores$dist_to_center),
            max(pca_scores$dist_to_center), dist_cv,
            n_outliers, outlier_ratio,
            median(rsd_values),
            sum(rsd_values < 20) / length(rsd_values) * 100,
            sum(rsd_values < 30) / length(rsd_values) * 100)
)
qc_report <- as.data.frame(summary_list)
write_csv(qc_report, "output/QC_assessment_summary.csv")
cat("📁 QC assessment summary saved: output/QC_assessment_summary.csv\n")

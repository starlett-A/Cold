# =================================================================
# 🌋 VOLCANO PLOT (PUBLICATION-QUALITY)
# Purpose: Generate optimised volcano plots for differential
#          metabolites, with intelligent label placement.
# Input:   Mapped significant metabolites from univariate analysis
# Output:  Multiple volcano plot variants (PNG & PDF)
# =================================================================

cat("🚀 Generating volcano plots...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("ggplot2", "ggrepel", "dplyr", "scales", "readr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD & PREPARE DATA
# =============================================
cat("\n📂 STEP 1: Loading significant metabolites (mapped)...\n")

input_file <- "output/univariate_significant_metabolites_Mapped.csv"
if (!file.exists(input_file)) {
  stop("Error: Mapped significant metabolites not found at: ", input_file,
       "\nPlease run 08_univariate_analysis.R then 10_name_mapping.R first.")
}

sig_df <- read_csv(input_file, show_col_types = FALSE)

# Check required columns
required_cols <- c("Metabolite", "p_adj_BH", "log2FC")
missing_cols <- setdiff(required_cols, colnames(sig_df))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

volcano_data <- sig_df %>%
  mutate(
    Metabolite = as.character(Metabolite),
    p_adj_BH   = as.numeric(p_adj_BH),
    log2FC     = as.numeric(log2FC),
    neg_log10_p = -log10(p_adj_BH),
    # Handle infinite values
    neg_log10_p = ifelse(is.infinite(neg_log10_p) | is.na(neg_log10_p),
                         max(neg_log10_p[!is.infinite(neg_log10_p)], na.rm = TRUE),
                         neg_log10_p),
    direction = ifelse(log2FC > 0, "Up in Poor", "Down in Poor")
  )

cat(sprintf("  Metabolites: %d\n", nrow(volcano_data)))
cat(sprintf("  Significant (FDR < 0.05): %d\n", sum(volcano_data$p_adj_BH < 0.05, na.rm = TRUE)))

# =============================================
# 2. VOLCANO PLOT FUNCTION
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 2: Creating optimised volcano function\n")
cat(paste0(strrep("=", 50), "\n"))

create_volcano <- function(data,
                           pval_threshold = 0.05,
                           fc_threshold   = 1.5,
                           top_n_labels   = 11,
                           title          = "Volcano plot: Fit-Good vs Fit-Poor",
                           point_size     = 5,
                           label_size     = 7) {
  
  plot_data <- data %>%
    mutate(
      significance = case_when(
        p_adj_BH < pval_threshold & log2FC > log2(fc_threshold) ~ "Up in Poor",
        p_adj_BH < pval_threshold & log2FC < -log2(fc_threshold) ~ "Up in Good",
        p_adj_BH < pval_threshold ~ "FDR only",
        abs(log2FC) > log2(fc_threshold) ~ "FC only",
        TRUE ~ "Not significant"
      ),
      label_priority = case_when(
        p_adj_BH < pval_threshold & abs(log2FC) > log2(fc_threshold) ~ "both",
        p_adj_BH < pval_threshold ~ "p",
        abs(log2FC) > log2(fc_threshold) ~ "fc",
        TRUE ~ "none"
      ),
      signif_score = -log10(p_adj_BH) * abs(log2FC)
    )
  
  # Select top N labels (prioritise both-sig, then p, then fc)
  labels_to_show <- plot_data %>%
    filter(label_priority != "none") %>%
    arrange(match(label_priority, c("both", "p", "fc")), -signif_score) %>%
    head(top_n_labels)
  
  cat(sprintf("  %d labels selected\n", nrow(labels_to_show)))
  
  # Build plot
  p <- ggplot(plot_data, aes(x = log2FC, y = neg_log10_p)) +
    # Non-significant points
    geom_point(data = filter(plot_data, significance == "Not significant"),
               color = "grey80", alpha = 0.3, size = point_size) +
    # Significant points
    geom_point(data = filter(plot_data, significance != "Not significant"),
               aes(color = significance), alpha = 1.0, size = point_size + 0.5) +
    # Threshold lines
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5, alpha = 0.5) +
    # Highlight labelled points
    geom_point(data = labels_to_show, shape = 21, color = "black", fill = NA,
               size = point_size + 2, stroke = 1.2) +
    # Labels
    geom_text_repel(
      data = labels_to_show,
      aes(label = Metabolite),
      size = label_size, fontface = "bold",
      max.overlaps = 100, min.segment.length = 0.3, segment.size = 0.6,
      segment.color = "grey40", segment.alpha = 0.7,
      box.padding = 0.8, point.padding = 0.6, force = 1.2, force_pull = 0.8,
      max.iter = 5000, nudge_x = 0.1, nudge_y = 0.1, direction = "both",
      seed = 12345, show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("Up in Poor" = "#E41A1C", "Up in Good" = "#377EB8",
                 "FDR only" = "#CCCCCC", "FC only" = "#984EA3"),
      name = "Significance"
    ) +
    labs(
      title = title,
      x = bquote(bold(log[2]~"(Fold Change, Poor/Good)")),
      y = bquote(bold(-log[10]~"(FDR-adjusted p-value)"))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 30, margin = margin(b = 20)),
      axis.title.x = element_text(face = "bold", size = 26, margin = margin(t = 15)),
      axis.title.y = element_text(face = "bold", size = 26, margin = margin(r = 15)),
      axis.text = element_text(size = 20, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.ticks.length = unit(0.25, "cm"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 24, margin = margin(b = 10)),
      legend.text = element_text(size = 22, color = "black"),
      legend.key.size = unit(1.2, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(1.2, "cm"),
      legend.spacing = unit(0.8, "cm"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
      plot.margin = margin(25, 25, 25, 25)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))
  
  return(list(plot = p, data = plot_data, labels = labels_to_show))
}

# =============================================
# 3. GENERATE MAIN VOLCANO PLOT
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 3: Creating main volcano plot\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

volcano_main <- create_volcano(
  data        = volcano_data,
  pval_threshold = 0.05,
  fc_threshold   = 1.5,
  top_n_labels   = 11,
  title          = "Volcano plot: Fit-Good vs Fit-Poor",
  point_size     = 4,
  label_size     = 5.5
)

ggsave("output/figures/volcano_optimised.png",
       volcano_main$plot, width = 16, height = 12, dpi = 600, bg = "white")
ggsave("output/figures/volcano_optimised.pdf",
       volcano_main$plot, width = 16, height = 12)
cat("✅ Main volcano plot saved\n")

# =============================================
# 4. ALTERNATIVE VERSIONS
# =============================================
cat("\n📊 Generating alternative volcano variants...\n")

# Variant 1: fewer labels (top 10)
volcano_10 <- create_volcano(volcano_data, top_n_labels = 10,
                             title = "Volcano plot (Top 10)",
                             point_size = 4, label_size = 5.5)
ggsave("output/figures/volcano_top10.png", volcano_10$plot,
       width = 16, height = 12, dpi = 600, bg = "white")
cat("✅ Top-10 version saved\n")

# Variant 2: smaller font
volcano_small <- create_volcano(volcano_data, top_n_labels = 11,
                                title = "Volcano plot (Smaller labels)",
                                point_size = 4, label_size = 4.5)
ggsave("output/figures/volcano_small_labels.png", volcano_small$plot,
       width = 16, height = 12, dpi = 600, bg = "white")
cat("✅ Small-label version saved\n")

# Variant 3: only both-significant labels
volcano_both <- create_volcano(volcano_data, top_n_labels = 11,
                               title = "Volcano plot (Both sig. only)")
# Restrict to only 'both' labels
volcano_both$labels <- volcano_both$labels %>%
  filter(p_adj_BH < 0.05 & abs(log2FC) > log2(1.5))
volcano_both_plot <- volcano_both$plot
ggsave("output/figures/volcano_both_significant.png",
       volcano_both_plot, width = 16, height = 12, dpi = 600, bg = "white")
cat("✅ Both-significant-only version saved\n")

# Save underlying data
write_csv(volcano_main$data, "output/volcano_plot_data.csv")
write_csv(volcano_main$labels, "output/volcano_plot_label_data.csv")
cat("✅ Volcano data files saved: output/volcano_plot_data.csv, output/volcano_plot_label_data.csv\n")

# =============================================
# 5. FINAL SUMMARY
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 VOLCANO PLOTS COMPLETE\n")
cat(paste0(strrep("=", 60), "\n\n"))

cat("📁 Generated files in output/figures/:\n")
cat("  - volcano_optimised.png / pdf\n")
cat("  - volcano_top10.png\n")
cat("  - volcano_small_labels.png\n")
cat("  - volcano_both_significant.png\n")
cat("Data exported to output/volcano_plot_data.csv\n")
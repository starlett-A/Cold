# =================================================================
# 🫧 PATHWAY ENRICHMENT BUBBLE PLOT
# Purpose: Visualise pathway enrichment results (significance & hits)
# Input:   Pre-defined pathway data (embedded)
# Output:  Publication-quality bubble plot (PNG & PDF)
# =================================================================

cat("🚀 Generating pathway enrichment bubble plot...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("ggplot2", "ggrepel"))
cat("✅ Package loading complete\n")

# =============================================
# 1. DEFINE PATHWAY DATA
# =============================================
cat("\n📊 STEP 1: Preparing pathway enrichment data...\n")

pathway_data <- data.frame(
  Pathway = c(
    "Tryptophan metabolism",
    "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "Nitrogen metabolism",
    "Purine metabolism",
    "Phenylalanine metabolism",
    "Arginine biosynthesis",
    "Butanoate metabolism",
    "Histidine metabolism",
    "Pantothenate and CoA biosynthesis"
  ),
  PValue = c(0.012525, 0.017489, 0.026134, 0.034686, 0.034714,
             0.060068, 0.064237, 0.068391, 0.084848),
  Hits   = c(2, 1, 1, 2, 1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)

# Categorise significance
pathway_data$Group <- ifelse(pathway_data$PValue < 0.05, "Significant", "Trend")

# Order pathways by P-value (descending for vertical display)
pathway_data$Pathway <- factor(
  pathway_data$Pathway,
  levels = pathway_data$Pathway[order(pathway_data$PValue)]
)

cat(sprintf("  Pathways: %d (Significant: %d, Trend: %d)\n",
            nrow(pathway_data),
            sum(pathway_data$Group == "Significant"),
            sum(pathway_data$Group == "Trend")))

# =============================================
# 2. CREATE BUBBLE PLOT
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 2: Building publication-quality bubble plot\n")
cat(paste0(strrep("=", 50), "\n"))

# Calculate x-axis range for -log10(P)
log_p <- -log10(pathway_data$PValue)
x_min <- floor(min(log_p) * 10) / 10
x_max <- ceiling(max(log_p) * 10) / 10

p <- ggplot(pathway_data, aes(x = -log10(PValue), y = Pathway)) +
  geom_point(aes(size = Hits, color = Group), alpha = 0.7) +
  
  # Significance threshold (P=0.05)
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "red", alpha = 0.5) +
  
  # P-value labels
  geom_text_repel(
    aes(label = paste0("P = ", format(PValue, digits = 3))),
    size = 5,
    box.padding = 0.3,
    color = "gray30",
    fontface = "bold"
  ) +
  
  # Colour & size scales
  scale_color_manual(
    values = c("Significant" = "#E41A1C", "Trend" = "#377EB8"),
    name = "Significance"
  ) +
  scale_size_continuous(
    range = c(4, 10),
    breaks = c(1, 2),
    name = "Number of\nMetabolites"
  ) +
  
  # X-axis breaks
  scale_x_continuous(
    limits = c(x_min - 0.1, x_max + 0.1),
    breaks = seq(x_min, x_max, by = 0.5),
    expand = c(0.02, 0.02)
  ) +
  
  labs(
    title = "Metabolic Pathway Enrichment",
    x = expression(bold(-log[10]*"(P-value)")),
    y = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5,
                              margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", size = 20, color = "black",
                                margin = margin(t = 10)),
    axis.text.x  = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y  = element_text(size = 15, color = "black"),
    
    legend.title = element_text(face = "bold", size = 20),
    legend.text  = element_text(size = 18),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.key.size = unit(1.8, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(1.2, "cm"),
    legend.spacing.y = unit(0.8, "cm"),
    
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  
  guides(color = guide_legend(override.aes = list(size = 6, shape = 19)))

# =============================================
# 3. SAVE FIGURES
# =============================================
cat("\n💾 Saving bubble plot...\n")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

ggsave("output/figures/pathway_enrichment_bubble.png",
       plot = p, width = 14, height = 10, dpi = 600, bg = "white")
ggsave("output/figures/pathway_enrichment_bubble.pdf",
       plot = p, width = 12, height = 10)

cat("✅ Plots saved in output/figures/\n")

cat(paste0("\n", strrep("=", 50), "\n"))
cat("🎉 PATHWAY BUBBLE PLOT COMPLETE\n")
cat(paste0(strrep("=", 50), "\n"))
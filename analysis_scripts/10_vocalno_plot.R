# =================================================================
# 🌋 VOLCANO PLOT (OPTIMISED, SPLIT LABELS)
# Purpose: Generate publication‑ready volcano plot with labels
#          pushed to left/right sides to avoid overlap.
# Input:   Mapped significant metabolites
# Output:  High‑resolution volcano plot (PNG & PDF)
# =================================================================

cat("🚀 Generating optimised volcano plot (split labels)...\n")
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
# 1. LOAD DATA
# =============================================
cat("\n📂 STEP 1: Loading mapped significant metabolites...\n")

input_file <- "output/univariate_significant_metabolites_Mapped.csv"
if (!file.exists(input_file)) {
  stop("Error: Mapped significant metabolites not found at: ", input_file,
       "\nPlease run 08_univariate_analysis.R then 10_name_mapping.R first.")
}

sig_data <- read_csv(input_file, show_col_types = FALSE)
cat(sprintf("  Metabolites loaded: %d\n", nrow(sig_data)))

# =============================================
# 2. DATA PREPROCESSING & LABEL SHORTENING
# =============================================
cat("\n🔧 STEP 2: Preparing data and creating short labels...\n")

volcano_data <- sig_data %>%
  mutate(
    Metabolite  = as.character(Metabolite),
    p_adj_BH    = as.numeric(p_adj_BH),
    log2FC      = as.numeric(log2FC),
    neg_log10_p = -log10(p_adj_BH),

    # Shorten metabolite names for cleaner labels
    short_label = case_when(
      Metabolite == "X4...2.4.Dihydroxy.3.3.dimethylbutanoyl.amino.butanoic.acid" ~ "Pantogab",
      Metabolite == "X4.11.13.15.Tetrahydroridentin.B" ~ "Tetrahydroridentin B",
      Metabolite == "N.Ornithyl.L.taurine" ~ "N-Ornithyl-L-taurine",
      Metabolite == "X9.Decenoylcarnitine" ~ "9-Decenoylcarnitine",
      Metabolite == "X8.Hydroxyadenine" ~ "8-Hydroxyadenine",
      Metabolite == "X2.Furoylglycine" ~ "2-Furoylglycine",
      Metabolite == "X11.Epicortisol" ~ "11-Epicortisol",
      Metabolite == "X12.Hydroxystearic.acid" ~ "12-Hydroxystearic acid",
      Metabolite == "piperazine.2.3.dione" ~ "Piperazine-2,3-dione",
      Metabolite == "Acetylenedicarboxylic.acid" ~ "Acetylenedicarboxylic acid",
      Metabolite == "Phenylalanine.betaine" ~ "Phenylalanine betaine",
      Metabolite == "Phenylalanylphenylalanine" ~ "Phe-Phe",
      Metabolite == "Pantothenic.acid" ~ "Pantothenic acid",
      Metabolite == "Glutamic.acid" ~ "Glutamic acid",
      Metabolite == "Indoleacetaldehyde" ~ "Indoleacetaldehyde",
      Metabolite == "o.Tyrosine" ~ "o-Tyrosine",
      Metabolite == "Phenylalanine" ~ "Phenylalanine",
      Metabolite == "Eucaglobulin" ~ "Eucaglobulin",
      Metabolite == "Methylhistidine" ~ "Methylhistidine",
      Metabolite == "Caprylic.acid" ~ "Caprylic acid",
      Metabolite == "Guanine" ~ "Guanine",
      Metabolite == "Serotonin" ~ "Serotonin",
      TRUE ~ Metabolite
    ),

    significance_score = neg_log10_p * abs(log2FC),

    significance = case_when(
      p_adj_BH < 0.05 & log2FC >  log2(1.5) ~ "Up in Fit-Poor",
      p_adj_BH < 0.05 & log2FC < -log2(1.5) ~ "Up in Fit-Good",
      p_adj_BH < 0.05                        ~ "FDR significant only",
      TRUE                                    ~ "Not significant"
    )
  )

# =============================================
# 3. SELECT METABOLITES TO LABEL
# =============================================
cat("\n🔖 STEP 3: Selecting metabolites for labelling (both FDR & FC significant)...\n")

labels_to_show <- volcano_data %>%
  filter(significance %in% c("Up in Fit-Poor", "Up in Fit-Good")) %>%
  arrange(desc(significance_score))

labels_left  <- labels_to_show %>% filter(log2FC < 0)   # Up in Fit-Good → left side
labels_right <- labels_to_show %>% filter(log2FC >= 0)  # Up in Fit-Poor → right side

cat(sprintf("  Total labelled: %d (left: %d, right: %d)\n",
            nrow(labels_to_show), nrow(labels_left), nrow(labels_right)))

# =============================================
# 4. BUILD VOLCANO PLOT (SPLIT‑SIDE LABELS)
# =============================================
cat("\n🎨 STEP 4: Creating volcano plot with split labels...\n")

p <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_p)) +

  # 4.1 Non‑significant points
  geom_point(data = filter(volcano_data, significance == "Not significant"),
             color = "grey80", alpha = 0.4, size = 4) +

  # 4.2 Significant points
  geom_point(data = filter(volcano_data, significance != "Not significant"),
             aes(color = significance), alpha = 1.0, size = 4.5) +

  # 4.3 Threshold lines
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4, alpha = 0.5) +

  # 4.4 Highlight labelled points
  geom_point(data = labels_to_show, shape = 21, color = "black", fill = NA,
             size = 6.5, stroke = 1.2) +

  # 4.5a Left labels (Up in Fit-Good) – push to left
  geom_text_repel(
    data          = labels_left,
    aes(label     = short_label),
    size          = 4.5,
    fontface      = "bold",
    color         = "#377EB8",
    direction     = "y",
    nudge_x       = -0.8,
    hjust         = 1,
    segment.size  = 0.5,
    segment.color = "grey40",
    segment.alpha = 0.7,
    box.padding   = 0.5,
    point.padding = 0.4,
    force         = 2,
    max.overlaps  = 50,
    max.iter      = 10000,
    seed          = 42,
    show.legend   = FALSE
  ) +

  # 4.5b Right labels (Up in Fit-Poor) – push to right
  geom_text_repel(
    data          = labels_right,
    aes(label     = short_label),
    size          = 4.5,
    fontface      = "bold",
    color         = "#E41A1C",
    direction     = "y",
    nudge_x       = 0.8,
    hjust         = 0,
    segment.size  = 0.5,
    segment.color = "grey40",
    segment.alpha = 0.7,
    box.padding   = 0.5,
    point.padding = 0.4,
    force         = 2,
    max.overlaps  = 50,
    max.iter      = 10000,
    seed          = 42,
    show.legend   = FALSE
  ) +

  # 4.6 Colour scale
  scale_color_manual(
    values = c(
      "Up in Fit-Poor"       = "#E41A1C",
      "Up in Fit-Good"       = "#377EB8",
      "FDR significant only" = "#999999"
    ),
    name = "Significance"
  ) +

  # 4.7 Axis labels and title
  labs(
    title = "Volcano Plot: Fit-Good vs Fit-Poor",
    x     = bquote(bold(log[2]~"(Fold Change, Poor/Good)")),
    y     = bquote(bold(-log[10]~"(FDR-adjusted p-value)"))
  ) +

  # 4.8 Theme
  theme_minimal(base_size = 14) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 28,
                                    margin = margin(b = 20)),
    axis.title.x     = element_text(face = "bold", size = 24, margin = margin(t = 12)),
    axis.title.y     = element_text(face = "bold", size = 24, margin = margin(r = 12)),
    axis.text        = element_text(size = 18, face = "bold", color = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    legend.title     = element_text(face = "bold", size = 20),
    legend.text      = element_text(size = 18),
    legend.position  = "right",
    legend.key.size  = unit(1.1, "cm"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 1.2),
    plot.margin      = margin(25, 40, 25, 40)   # extra side margins for labels
  ) +

  # 4.9 Axis expansions (extra room for labels)
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.20)))  # 20% on each side

# =============================================
# 5. SAVE PLOTS
# =============================================
cat("\n💾 STEP 5: Saving volcano plot...\n")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

ggsave("output/figures/volcano_plot_final.png",
       p, width = 18, height = 12, dpi = 600, bg = "white")
ggsave("output/figures/volcano_plot_final.pdf",
       p, width = 18, height = 12)
cat("✅ Saved: output/figures/volcano_plot_final.png / .pdf\n")

cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 VOLCANO PLOT COMPLETE\n")
cat(paste0(strrep("=", 60), "\n"))

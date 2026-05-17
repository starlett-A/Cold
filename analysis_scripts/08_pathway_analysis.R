# =================================================================
# 🧬 METABOANALYST PATHWAY ENRICHMENT (ORA)
# Purpose: Over‑representation analysis (ORA) using MetaboAnalystR,
#          given a list of HMDB IDs from significant metabolites.
# Input:   HMDB IDs (embedded), or read from a file if preferred
# Output:  Enrichment table, bubble plot
# =================================================================

cat("🚀 Starting MetaboAnalystR pathway enrichment...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("MetaboAnalystR", "ggplot2", "ggrepel", "dplyr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. DEFINE HMDB IDs (from your significant metabolites)
# =============================================
cat("\n📂 STEP 1: Defining HMDB ID list...\n")

my_hmdb_ids <- c(
  "HMDB0000148",
  "HMDB0000439",
  "HMDB0000210",
  "HMDB0013302",
  "HMDB0247933",
  "HMDB0000542",
  "HMDB0256573",
  "HMDB0001190",
  "HMDB0006050",
  "HMDB0000259",
  "HMDB0061706",
  "HMDB0036742",
  "HMDB0000159",
  "HMDB0000479",
  "HMDB0246316",
  "HMDB0240552",
  "HMDB0036150",
  "HMDB0033519",
  "HMDB0000482",
  "HMDB0000132",
  "HMDB0013205",
  "HMDB0244354"
)

cat(sprintf("  Total HMDB IDs: %d\n", length(my_hmdb_ids)))

# Optional: read from an external file if you prefer
# my_hmdb_ids <- readLines("output/sig_hmdb_list.txt")

# =============================================
# 2. SET UP METABOANALYST OBJECT & MAP IDs
# =============================================
cat("\n🔗 STEP 2: Mapping HMDB IDs to KEGG...\n")

mSet <- InitDataObjects("conc", "msetqea", FALSE)
mSet <- Setup.MapData(mSet, my_hmdb_ids)
mSet <- CrossReferencing(mSet, "hmdb")
mSet <- CreateMappingResultTable(mSet)

# Show mapping summary
map_table <- mSet$dataSet$map.table
cat(sprintf("  Mapped metabolites: %d\n", nrow(map_table)))
print(head(map_table))

# =============================================
# 3. OVER‑REPRESENTATION ANALYSIS
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 3: Running ORA (hypergeometric test)\n")
cat(paste0(strrep("=", 50), "\n"))

# Use all KEGG metabolites as background
mSet <- SetMetabolomeFilter(mSet, F)

# Calculate ORA scores (method: "rbc" or "hyperg")
mSet <- CalculateOraScore(mSet, "rbc", "hyperg")

# Extract enrichment results
enrich_res <- mSet$analSet$ora.mat
if (is.null(enrich_res) || nrow(enrich_res) == 0) {
  stop("No enrichment results found. Please check your HMDB IDs.")
}
cat(sprintf("  Enriched pathways: %d\n", nrow(enrich_res)))

# Save full enrichment table
write.csv(enrich_res, "output/pathway_enrichment_results.csv", row.names = FALSE)
cat("✅ Full enrichment table saved: output/pathway_enrichment_results.csv\n")

# =============================================
# 4. PREPARE DATA FOR BUBBLE PLOT
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 4: Generating bubble plot\n")
cat(paste0(strrep("=", 50), "\n"))

# Extract pathway names, raw p-values, and hits
res <- data.frame(
  Pathway = rownames(enrich_res),
  PValue  = as.numeric(enrich_res[, "Raw p"]),
  Hits    = as.numeric(enrich_res[, "Hits"]),
  stringsAsFactors = FALSE
)

# Keep only pathways with P < 0.1
res <- res[res$PValue < 0.1 & complete.cases(res), ]
if (nrow(res) == 0) {
  cat("⚠️  No pathways with P < 0.1. Consider relaxing the threshold.\n")
  q(save = "no")
}
cat(sprintf("  Pathways with P < 0.1: %d\n", nrow(res)))

# Categorise significance
res$Group <- ifelse(res$PValue < 0.05, "Significant", "Trend")
res <- res[order(res$PValue), ]
res$Pathway <- factor(res$Pathway, levels = res$Pathway)

# X-axis range for -log10(P)
log_p <- -log10(res$PValue)
min_log_p <- floor(min(log_p) * 10) / 10
max_log_p <- ceiling(max(log_p) * 10) / 10

# Optional: map KEGG codes to human-readable names if desired
# You can customise this list as needed
kegg_to_name <- c(
  "hsa00340" = "Histidine metabolism",
  "hsa00400" = "Phenylalanine, tyrosine and tryptophan biosynthesis",
  "hsa00380" = "Tryptophan metabolism",
  "hsa00910" = "Nitrogen metabolism",
  "hsa00360" = "Phenylalanine metabolism",
  "hsa00230" = "Purine metabolism",
  "hsa00650" = "Butanoate metabolism",
  "hsa00220" = "Arginine biosynthesis"
)

# If the pathway codes match the KEGG map, replace them with readable names
# (otherwise keep original)
res$Display_Name <- ifelse(res$Pathway %in% names(kegg_to_name),
                           kegg_to_name[res$Pathway], as.character(res$Pathway))
res$Pathway <- factor(res$Display_Name, levels = unique(res$Display_Name))

# =============================================
# 5. CREATE BUBBLE PLOT
# =============================================
bubble_plot <- ggplot(res, aes(x = -log10(PValue), y = Pathway,
                               size = Hits, color = Group)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "red", alpha = 0.5) +
  geom_text_repel(aes(label = paste0("P=", format(PValue, digits = 3))),
                  size = 5, box.padding = 0.3, color = "gray30",
                  fontface = "bold") +
  scale_color_manual(values = c("Significant" = "#E41A1C",
                                "Trend" = "#377EB8")) +
  scale_size_continuous(range = c(4, 10),
                        breaks = sort(unique(res$Hits))) +
  scale_x_continuous(limits = c(min_log_p - 0.1, max_log_p + 0.1),
                     breaks = seq(min_log_p, max_log_p, by = 0.5)) +
  labs(title = "Metabolic Pathway Enrichment",
       x = expression(bold(-log[10] * "(P-value)")),
       y = NULL,
       size = "Number of\nMetabolites",
       color = "Significance") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5,
                              margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", size = 20, color = "black",
                                margin = margin(t = 10)),
    axis.text.x = element_text(size = 18, color = "black", face = "bold"),
    axis.text.y = element_text(size = 15, color = "black"),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18),
    legend.position = "right",
    legend.key.size = unit(1.8, "cm"),
    legend.spacing.y = unit(0.8, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 6, shape = 19)))

# =============================================
# 6. SAVE PLOT
# =============================================
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

ggsave("output/figures/metaboanalyst_pathway_bubble.png",
       plot = bubble_plot, width = 14, height = 10, dpi = 600)
ggsave("output/figures/metaboanalyst_pathway_bubble.pdf",
       plot = bubble_plot, width = 14, height = 10)
cat("✅ Bubble plot saved in output/figures/\n")

cat(paste0("\n", strrep("=", 50), "\n"))
cat("🎉 METABOANALYST PATHWAY ENRICHMENT COMPLETE\n")
cat(paste0(strrep("=", 50), "\n"))

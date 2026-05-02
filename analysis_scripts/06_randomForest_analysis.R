# =================================================================
# 🌲 RANDOM FOREST FEATURE SELECTION
# Purpose: Validate important metabolites with Random Forest
# Input:   Scaled data, VIP important list, Bootstrap results (optional)
# Output:  Importance scores, selected metabolites abundance,
#          integrated RF-Bootstrap table, diagnostic barplot
# =================================================================

cat("🚀 Starting Random Forest analysis...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("randomForest", "dplyr", "ggplot2", "readr", "magrittr", "tidyr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD INPUT DATA
# =============================================
cat("\n📂 STEP 1: Loading input files...\n")

# Scaled data (from PCA step)
scaled_file <- "output/preprocessed_data_scaled.csv"
if (!file.exists(scaled_file)) {
  stop("Error: Scaled data not found at: ", scaled_file,
       "\nPlease run 03_pca_analysis.R first.")
}
scaled_data <- read_csv(scaled_file, show_col_types = FALSE)
cat(sprintf("  Scaled data: %d samples, %d metabolites\n",
            nrow(scaled_data), ncol(scaled_data) - 2))

# VIP important list
vip_file <- "output/VIP_important.csv"
if (!file.exists(vip_file)) {
  stop("Error: VIP important table not found at: ", vip_file,
       "\nPlease run 04_plsda_analysis.R first.")
}
vip_df <- read_csv(vip_file, show_col_types = FALSE)
vip_features <- vip_df$Metabolite
cat(sprintf("  VIP important metabolites: %d\n", length(vip_features)))

# Name mapping (for figures)
mapping_file <- "output/metabolite_name_mapping.csv"
name_mapping <- NULL
if (file.exists(mapping_file)) {
  name_mapping <- read_csv(mapping_file, show_col_types = FALSE)
  cat("  Name mapping loaded\n")
} else {
  cat("  ⚠️  Name mapping not found; will use safe names in plots\n")
}

# =============================================
# 2. PREPARE DATA MATRIX
# =============================================
cat("\n🔧 STEP 2: Preparing feature matrix (VIP metabolites only)...\n")

# Ensure Sample_ID exists
if (!"Sample_ID" %in% colnames(scaled_data)) {
  stop("Scaled data must contain a 'Sample_ID' column")
}
samples <- scaled_data$Sample_ID
Groups  <- factor(scaled_data$Group, levels = c("Fit-Good", "Fit-Poor"))

# Keep only columns that exist and are in VIP list
available_features <- intersect(vip_features, colnames(scaled_data))
if (length(available_features) == 0) {
  stop("No VIP metabolites found in scaled data. Check metabolite names.")
}
cat(sprintf("  Usable VIP features: %d / %d\n",
            length(available_features), length(vip_features)))

X_vip <- as.data.frame(scaled_data[, available_features, drop = FALSE])
Y     <- Groups

# =============================================
# 3. RANDOM FOREST MODEL
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 3: Training Random Forest (ntree=1000)\n")
cat(paste0(strrep("=", 50), "\n"))

set.seed(2026)
rf_model <- randomForest(
  x           = X_vip,
  y           = Y,
  ntree       = 1000,
  importance  = TRUE,
  proximity   = TRUE,
  do.trace    = 50
)

oob_error <- mean(rf_model$err.rate[, "OOB"])
oob_accuracy <- (1 - oob_error) * 100
cat(sprintf("\n  OOB error rate: %.2f%%\n", oob_error * 100))
cat(sprintf("  OOB accuracy  : %.2f%%\n", oob_accuracy))

# =============================================
# 4. FEATURE IMPORTANCE
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 4: Extracting feature importance\n")
cat(paste0(strrep("=", 50), "\n"))

imp_vec <- rf_model$importance[, "MeanDecreaseAccuracy"]
imp_df <- data.frame(
  Metabolite     = names(imp_vec),
  RF_Importance  = imp_vec,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(RF_Importance)) %>%
  mutate(Rank = row_number())

# Select top 30% as RF-confirmed important
n_select <- max(1, round(0.3 * nrow(imp_df)))
rf_selected <- imp_df$Metabolite[1:n_select]

cat(sprintf("  Top 30%% selected: %d metabolites\n", n_select))

# Save importance table
write_csv(imp_df, "output/RandomForest_importance.csv")
cat("✅ Saved: output/RandomForest_importance.csv\n")

# =============================================
# 5. SAVE ABUNDANCE TABLE (RF selected)
# =============================================
cat("\n💾 Saving abundance data for RF-selected metabolites...\n")

abundance_df <- data.frame(
  Sample_ID = samples,
  Group     = as.character(Y),
  scaled_data[, rf_selected, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_csv(abundance_df, "output/RF_selected_metabolites_abundance.csv")
cat("✅ Saved: output/RF_selected_metabolites_abundance.csv\n")

# =============================================
# 6. INTEGRATE WITH BOOTSTRAP RESULTS (if available)
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 5: Integrating with Bootstrap stability...\n")
cat(paste0(strrep("=", 50), "\n"))

bootstrap_file <- "output/bootstrap_stability_results.csv"
if (file.exists(bootstrap_file)) {
  boot_df <- read_csv(bootstrap_file, show_col_types = FALSE)
  
  integrated <- imp_df %>%
    left_join(boot_df %>% select(Metabolite, Bootstrap_Frequency, Composite_Score, Confidence_Level),
              by = "Metabolite")
  
  # Normalise RF importance for combined score (if Bootstrap frequency is missing, skip)
  if (!all(is.na(integrated$Bootstrap_Frequency))) {
    scale01 <- function(x) {
      if (length(unique(na.omit(x))) <= 1) return(rep(0, length(x)))
      (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    }
    
    integrated <- integrated %>%
      mutate(
        RF_Normalized       = scale01(RF_Importance),
        Bootstrap_Normalized = ifelse(is.na(Bootstrap_Frequency), 0, Bootstrap_Frequency / 100),
        Combined_Score       = 0.6 * RF_Normalized + 0.4 * Bootstrap_Normalized,
        Combined_Rank        = rank(-Combined_Score, ties.method = "min")
      ) %>%
      arrange(desc(Combined_Score))
    
    write_csv(integrated, "output/Integrated_RF_Bootstrap.csv")
    cat("✅ Saved: output/Integrated_RF_Bootstrap.csv\n")
    
    # Show top 10
    cat("\n🏆 Top 10 by combined score:\n")
    top10 <- head(integrated, 10)
    for (i in 1:nrow(top10)) {
      cat(sprintf("  %2d. %-30s  RF=%.3f  Freq=%.1f%%  Combined=%.3f\n",
                  i, top10$Metabolite[i], top10$RF_Importance[i],
                  top10$Bootstrap_Frequency[i], top10$Combined_Score[i]))
    }
  } else {
    cat("⚠️  Bootstrap results missing Bootstrap_Frequency column. Integration skipped.\n")
  }
} else {
  cat("⚠️  Bootstrap results not found. Skipping integration.\n")
}

# =============================================
# 7. VISUALIZATION
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 6: Generating importance barplot\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# Prepare top 20 with display names (original if available)
plot_data <- imp_df %>%
  head(20)

if (!is.null(name_mapping)) {
  # Try to join on Safe_Name
  if ("Safe_Name" %in% colnames(name_mapping)) {
    plot_data <- plot_data %>%
      left_join(name_mapping, by = c("Metabolite" = "Safe_Name"))
  } else {
    # Fallback: assume first two columns are safe, original
    name_mapping <- name_mapping %>% rename(Safe_Name = 1, Original_Name = 2)
    plot_data <- plot_data %>%
      left_join(name_mapping, by = c("Metabolite" = "Safe_Name"))
  }
  plot_data <- plot_data %>%
    mutate(Display_Name = ifelse(is.na(Original_Name), Metabolite, Original_Name))
} else {
  plot_data <- plot_data %>% mutate(Display_Name = Metabolite)
}

# Order by importance
plot_data$Display_Name <- factor(
  plot_data$Display_Name,
  levels = plot_data$Display_Name[order(plot_data$RF_Importance)]
)

mean_imp <- mean(imp_df$RF_Importance, na.rm = TRUE)

p <- ggplot(plot_data, aes(x = Display_Name, y = RF_Importance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", RF_Importance)),
            hjust = -0.1, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = mean_imp, linetype = "dashed", color = "red", linewidth = 0.7) +
  coord_flip() +
  labs(
    title = "Top 20 Metabolites by Random Forest Importance",
    subtitle = sprintf("OOB Accuracy: %.1f%% | Top 30%% selected: %d",
                       oob_accuracy, n_select),
    x = "Metabolite",
    y = "Importance (Mean Decrease Accuracy)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.title         = element_text(size = 13, face = "bold"),
    plot.title         = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle      = element_text(color = "gray40", hjust = 0.5, size = 12),
    panel.grid.major   = element_line(color = "gray90"),
    panel.grid.minor   = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("output/figures/RandomForest_importance_top20.png", p, width = 12, height = 8, dpi = 600, bg = "white")
ggsave("output/figures/RandomForest_importance_top20.pdf", p, width = 12, height = 8)
cat("✅ Barplot saved\n")

# =============================================
# 8. FINAL SUMMARY
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 RANDOM FOREST ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 60), "\n\n"))

cat(sprintf("  OOB accuracy: %.1f%%\n", oob_accuracy))
cat(sprintf("  Top 30%% selected: %d metabolites\n", n_select))
cat("\n📁 Output files:\n")
cat("  output/RandomForest_importance.csv\n")
cat("  output/RF_selected_metabolites_abundance.csv\n")
if (file.exists("output/Integrated_RF_Bootstrap.csv")) {
  cat("  output/Integrated_RF_Bootstrap.csv\n")
}
cat("  output/figures/RandomForest_importance_top20.png/pdf\n")

cat("\n🔜 Next: combine results or proceed to pathway enrichment.\n")
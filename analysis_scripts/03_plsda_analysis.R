# =================================================================
# 🔬 PLS-DA ANALYSIS (SUPERVISED)
# Purpose: PLS-DA with LOOCV, permutation test, VIP selection
# Input:  Pareto-scaled data from 03_pca_analysis.R
# Output: Model results, VIP tables, performance metrics, figures
# =================================================================

cat("🚀 Starting PLS-DA analysis...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("mixOmics", "dplyr", "ggplot2", "tidyr", "readr", "magrittr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD SCALED DATA
# =============================================
cat("\n📂 STEP 1: Loading scaled data...\n")

# Input file from PCA step
input_file <- "output/preprocessed_data_scaled.csv"

if (!file.exists(input_file)) {
  stop("Error: Scaled data not found at: ", input_file,
       "\nPlease run 03_pca_analysis.R first (Pareto scaling).")
}

scaled_data <- read_csv(input_file, show_col_types = FALSE)

cat(sprintf("  Data loaded: %d samples, %d metabolites\n",
            nrow(scaled_data), ncol(scaled_data) - 2))
cat(sprintf("  Groups: %s\n", paste(unique(scaled_data$Group), collapse = ", ")))

# =============================================
# 2. PREPARE DATA MATRIX
# =============================================
cat("\n🔧 STEP 2: Preparing PLS-DA input...\n")

plsda_X <- as.matrix(scaled_data[, 3:ncol(scaled_data)])
rownames(plsda_X) <- scaled_data$Sample_ID
plsda_Y <- factor(scaled_data$Group, levels = c("Fit-Good", "Fit-Poor"))

cat("    Group distribution:\n")
print(table(plsda_Y))

# =============================================
# 3. MANUAL LOOCV FOR COMPONENT SELECTION
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("STEP 3: Leave-One-Out Cross-Validation\n")
cat(paste0(strrep("=", 60), "\n"))

# Test 1 to 3 components
error_rates <- numeric(3)
names(error_rates) <- paste0("comp", 1:3)

for (ncomp in 1:3) {
  cat(sprintf("  Testing %d component(s)... ", ncomp))
  predictions <- character(length(plsda_Y))
  
  for (i in seq_along(plsda_Y)) {
    train_X <- plsda_X[-i, , drop = FALSE]
    train_Y <- plsda_Y[-i]
    test_X  <- plsda_X[i, , drop = FALSE]
    
    if (length(unique(train_Y)) < 2) {
      predictions[i] <- as.character(train_Y[1])
    } else {
      model <- plsda(train_X, train_Y, ncomp = ncomp)
      pred  <- predict(model, test_X)$class$max.dist[, ncomp]
      predictions[i] <- as.character(pred)
    }
  }
  
  error_rates[ncomp] <- 1 - mean(predictions == as.character(plsda_Y))
  cat(sprintf("error = %.3f (accuracy = %.1f%%)\n",
              error_rates[ncomp], (1 - error_rates[ncomp]) * 100))
}

# Select optimal ncomp (lowest error, smallest ncomp in ties)
optimal_ncomp <- which.min(error_rates)
if (length(optimal_ncomp) > 1) optimal_ncomp <- min(optimal_ncomp)

cat(paste0("\n", strrep("=", 60), "\n"))
cat("LOOCV Summary\n")
cat(paste0(strrep("=", 60), "\n"))
for (i in seq_along(error_rates)) {
  star <- ifelse(i == optimal_ncomp, " ★ (selected)", "")
  cat(sprintf("  ncomp = %d: error = %.3f, accuracy = %.1f%%%s\n",
              i, error_rates[i], (1 - error_rates[i]) * 100, star))
}
cat(sprintf("\n✅ Optimal number of components: %d\n", optimal_ncomp))

# =============================================
# 4. FINAL PLS-DA MODEL
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 4: Fitting final PLS-DA model\n")
cat(paste0(strrep("=", 50), "\n"))

plsda_result <- mixOmics::plsda(plsda_X, plsda_Y, ncomp = optimal_ncomp)
cat("✅ Final model fitted\n")

# =============================================
# 5. VIP CALCULATION
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 5: Calculating VIP scores\n")
cat(paste0(strrep("=", 50), "\n"))

vip_matrix <- mixOmics::vip(plsda_result)

if (optimal_ncomp > 1) {
  avg_vip <- apply(vip_matrix[, 1:optimal_ncomp], 1, mean)
} else {
  avg_vip <- vip_matrix[, 1]
}

vip_df <- data.frame(
  Metabolite = colnames(plsda_X),
  VIP_Score  = round(avg_vip, 4),
  stringsAsFactors = FALSE
) %>% arrange(desc(VIP_Score))

vip_threshold <- 1.5
important_metabolites <- vip_df %>% filter(VIP_Score > vip_threshold)

cat(sprintf("  VIP > %.1f : %d metabolites\n", vip_threshold, nrow(important_metabolites)))
cat(sprintf("  VIP > 2.0 : %d metabolites\n", sum(vip_df$VIP_Score > 2.0)))
cat(sprintf("  VIP > 2.5 : %d metabolites\n", sum(vip_df$VIP_Score > 2.5)))

# =============================================
# 6. PERMUTATION TEST
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 6: Permutation test (n = 1000)\n")
cat(paste0(strrep("=", 50), "\n"))

perform_robust_permutation_test <- function(X, Y, ncomp, n_perm = 1000, seed = 2026) {
  set.seed(seed)
  n_samples <- nrow(X)
  
  # Original model LOOCV accuracy
  original_predictions <- character(n_samples)
  pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
  for (i in 1:n_samples) {
    train_idx <- setdiff(1:n_samples, i)
    model <- mixOmics::plsda(X[train_idx, , drop = FALSE], Y[train_idx], ncomp = ncomp)
    pred  <- predict(model, X[i, , drop = FALSE])$class$max.dist[, ncomp]
    original_predictions[i] <- as.character(pred)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  original_accuracy <- mean(original_predictions == as.character(Y))
  
  # Permuted labels
  perm_accuracies <- numeric(n_perm)
  pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
  for (p in 1:n_perm) {
    Y_perm <- sample(Y)
    perm_predictions <- character(n_samples)
    
    for (i in 1:n_samples) {
      train_idx <- setdiff(1:n_samples, i)
      if (length(unique(Y_perm[train_idx])) < 2) {
        perm_predictions[i] <- NA
        next
      }
      model_perm <- mixOmics::plsda(X[train_idx, , drop = FALSE], Y_perm[train_idx], ncomp = ncomp)
      pred_perm  <- predict(model_perm, X[i, , drop = FALSE])$class$max.dist[, ncomp]
      perm_predictions[i] <- as.character(pred_perm)
    }
    
    valid <- !is.na(perm_predictions)
    if (sum(valid) > 0) {
      perm_accuracies[p] <- mean(perm_predictions[valid] == as.character(Y_perm)[valid], na.rm = TRUE)
    } else {
      perm_accuracies[p] <- 0.5
    }
    setTxtProgressBar(pb, p)
  }
  close(pb)
  
  p_value <- sum(perm_accuracies >= original_accuracy) / n_perm
  
  cat(sprintf("\n  Original accuracy: %.4f (%.1f%%)\n", original_accuracy, original_accuracy * 100))
  cat(sprintf("  Mean permuted accuracy: %.4f\n", mean(perm_accuracies, na.rm = TRUE)))
  cat(sprintf("  Permutation p-value: %.4f\n", p_value))
  
  list(original_accuracy = original_accuracy,
       perm_accuracies   = perm_accuracies,
       p_value           = p_value,
       n_perm            = n_perm)
}

perm_results <- perform_robust_permutation_test(
  X = plsda_X,
  Y = plsda_Y,
  ncomp = optimal_ncomp,
  n_perm = 1000,
  seed = 2026
)

# =============================================
# 7. SAVE RESULTS
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 7: Saving results tables\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

# 7.1 Permutation results
perm_acc_df <- data.frame(
  Permutation_ID = seq_along(perm_results$perm_accuracies),
  Accuracy = perm_results$perm_accuracies
)
write_csv(perm_acc_df, "output/permutation_accuracies.csv")

perm_summary <- data.frame(
  Metric = c("Original_Accuracy", "P_Value", "N_Permutations",
             "Mean_Perm_Accuracy", "SD_Perm_Accuracy", "Median_Perm_Accuracy",
             "Significant_0.05", "Significant_0.01"),
  Value  = c(perm_results$original_accuracy,
             perm_results$p_value,
             perm_results$n_perm,
             mean(perm_results$perm_accuracies),
             sd(perm_results$perm_accuracies),
             median(perm_results$perm_accuracies),
             perm_results$p_value < 0.05,
             perm_results$p_value < 0.01)
)
write_csv(perm_summary, "output/permutation_summary.csv")
cat("✅ Permutation results saved: output/permutation_accuracies.csv, output/permutation_summary.csv\n")

# 7.2 VIP tables
write_csv(vip_df, "output/VIP_all.csv")
write_csv(important_metabolites, "output/VIP_important.csv")
write_csv(head(vip_df, 50), "output/VIP_top50.csv")
cat("✅ VIP tables saved: output/VIP_all.csv, output/VIP_important.csv, output/VIP_top50.csv\n")

# 7.3 Performance metrics
performance_table <- data.frame(
  Category = c("Model", "Model", "Validation", "Feature Selection",
               "Feature Selection", "Feature Selection", "Feature Selection"),
  Metric   = c("Optimal components", "LOOCV accuracy", "Permutation p-value",
               "VIP threshold", "VIP > 1.5", "VIP > 2.0", "VIP > 2.5"),
  Value    = c(optimal_ncomp,
               sprintf("%.1f%%", (1 - error_rates[optimal_ncomp]) * 100),
               sprintf("%.4f", perm_results$p_value),
               "1.5",
               sprintf("%d (%.1f%%)", nrow(important_metabolites),
                       nrow(important_metabolites) / nrow(vip_df) * 100),
               sprintf("%d (%.1f%%)", sum(vip_df$VIP_Score > 2.0),
                       sum(vip_df$VIP_Score > 2.0) / nrow(vip_df) * 100),
               sprintf("%d (%.1f%%)", sum(vip_df$VIP_Score > 2.5),
                       sum(vip_df$VIP_Score > 2.5) / nrow(vip_df) * 100))
)
write_csv(performance_table, "output/PLSDA_performance_metrics.csv")
cat("✅ Performance metrics saved: output/PLSDA_performance_metrics.csv\n")

# 7.4 Component 1 scores (if ncomp = 1)
if (optimal_ncomp == 1) {
  sample_scores <- data.frame(
    Sample_ID = rownames(plsda_X),
    Group = plsda_Y,
    Comp1_Score = plsda_result$variates$X[, 1]
  )
  write_csv(sample_scores, "output/PLSDA_component1_scores.csv")
  cat("✅ Component 1 scores saved: output/PLSDA_component1_scores.csv\n")
}

# =============================================
# 8. GENERATE FIGURES
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 8: Generating publication figures\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# 8.1 VIP distribution histogram
vip_dist_plot <- ggplot(vip_df, aes(x = VIP_Score)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = c(1.5, 2.0, 2.5),
             linetype = "dashed", color = c("red", "orange", "darkgreen"),
             linewidth = 0.8) +
  annotate("text", x = c(1.5, 2.0, 2.5), y = Inf,
           label = c("VIP>1.5", "VIP>2.0", "VIP>2.5"),
           color = c("red", "orange", "darkgreen"),
           hjust = -0.1, vjust = 1.5, size = 3.5) +
  labs(title = "Distribution of VIP Scores",
       subtitle = sprintf("Optimal ncomp = %d | Total metabolites: %d",
                          optimal_ncomp, nrow(vip_df)),
       x = "VIP Score", y = "Number of metabolites") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        panel.grid.minor = element_blank())

ggsave("output/figures/VIP_distribution.png", vip_dist_plot, width = 10, height = 6, dpi = 600)
ggsave("output/figures/VIP_distribution.pdf", vip_dist_plot, width = 10, height = 6)
cat("✅ VIP distribution plot saved\n")

# 8.2 Component selection plot
comp_sel_df <- data.frame(
  ncomp = 1:length(error_rates),
  Accuracy = (1 - error_rates) * 100
)

comp_sel_plot <- ggplot(comp_sel_df, aes(x = ncomp, y = Accuracy)) +
  geom_line(color = "grey60", linewidth = 1.2) +
  geom_point(aes(color = ifelse(ncomp == optimal_ncomp, "Selected", "Other")),
             size = 8, shape = 18) +
  geom_point(size = 4, color = "white") +
  scale_color_manual(values = c("Selected" = "red", "Other" = "steelblue"),
                     guide = "none") +
  scale_x_continuous(breaks = 1:length(error_rates)) +
  labs(title = "Determination of Optimal PLS-DA Components",
       subtitle = "Leave-One-Out Cross-Validation (LOOCV)",
       x = "Number of latent components",
       y = "LOOCV Accuracy (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()) +
  annotate("text", x = optimal_ncomp,
           y = comp_sel_df$Accuracy[optimal_ncomp] + 3,
           label = sprintf("Selected\nncomp = %d\nAccuracy = %.1f%%",
                           optimal_ncomp, comp_sel_df$Accuracy[optimal_ncomp]),
           color = "red", size = 3.5)

ggsave("output/figures/PLSDA_component_selection.png", comp_sel_plot, width = 8, height = 8, dpi = 600)
ggsave("output/figures/PLSDA_component_selection.pdf", comp_sel_plot, width = 8, height = 8)
cat("✅ Component selection plot saved\n")

# 8.3 Permutation histogram
perm_hist <- ggplot(data.frame(Accuracy = perm_results$perm_accuracies), aes(x = Accuracy)) +
  geom_histogram(bins = 30, fill = "lightblue", alpha = 0.8, color = "black") +
  geom_vline(xintercept = perm_results$original_accuracy,
             color = "red", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = perm_results$original_accuracy, y = Inf,
           label = sprintf("Original\nAccuracy = %.1f%%", perm_results$original_accuracy * 100),
           color = "red", hjust = -0.1, vjust = 1.5, size = 4) +
  labs(title = "Permutation Test Results",
       subtitle = sprintf("%d permutations, p = %.4f", perm_results$n_perm, perm_results$p_value),
       x = "Accuracy of models with permuted labels",
       y = "Frequency") +
  theme_minimal()

ggsave("output/figures/Permutation_histogram.png", perm_hist, width = 10, height = 6, dpi = 600)
ggsave("output/figures/Permutation_histogram.pdf", perm_hist, width = 10, height = 6)
cat("✅ Permutation histogram saved\n")

# 8.4 Top 20 VIP barplot
top_20 <- head(vip_df, 20)
vip_bar <- ggplot(top_20, aes(x = reorder(Metabolite, VIP_Score), y = VIP_Score)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Top 20 Metabolites by VIP Score",
       subtitle = sprintf("Optimal ncomp = %d, Accuracy = %.1f%%",
                          optimal_ncomp, (1 - error_rates[optimal_ncomp]) * 100),
       x = "Metabolite", y = "VIP Score") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))

ggsave("output/figures/Top20_VIP.png", vip_bar, width = 12, height = 8, dpi = 600)
ggsave("output/figures/Top20_VIP.pdf", vip_bar, width = 12, height = 8)
cat("✅ Top 20 VIP barplot saved\n")

# Special plot for ncomp = 1: Component 1 boxplot
if (optimal_ncomp == 1) {
  sample_scores <- read_csv("output/PLSDA_component1_scores.csv", show_col_types = FALSE)
  
  comp1_box <- ggplot(sample_scores, aes(x = Group, y = Comp1_Score, fill = Group)) +
    geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = Group), width = 0.15, size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 5, fill = "yellow", color = "black") +
    scale_fill_manual(values = c("Fit-Good" = "#2E86AB", "Fit-Poor" = "#A23B72")) +
    scale_color_manual(values = c("Fit-Good" = "#2E86AB", "Fit-Poor" = "#A23B72"), guide = "none") +
    labs(title = "Group Separation on PLS-DA Component 1",
         subtitle = paste0("ncomp = ", optimal_ncomp),
         x = NULL, y = "Component 1 Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # Add t-test p-value
  t_test <- t.test(Comp1_Score ~ Group, data = sample_scores)
  p_text <- ifelse(t_test$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", t_test$p.value))
  y_max <- max(sample_scores$Comp1_Score) * 1.15
  comp1_box <- comp1_box +
    annotate("segment", x = 1, xend = 2, y = y_max, yend = y_max, linewidth = 0.8) +
    annotate("text", x = 1.5, y = y_max * 1.02, label = p_text, size = 5, fontface = "bold")
  
  ggsave("output/figures/PLSDA_Comp1_boxplot.png", comp1_box, width = 7, height = 6, dpi = 600)
  ggsave("output/figures/PLSDA_Comp1_boxplot.pdf", comp1_box, width = 7, height = 6)
  cat("✅ Component 1 boxplot saved (ncomp=1)\n")
}

# =============================================
# 9. FINAL SUMMARY
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 PLS-DA ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 60), "\n\n"))

cat(sprintf("  Optimal components : %d\n", optimal_ncomp))
cat(sprintf("  LOOCV accuracy      : %.1f%%\n", (1 - error_rates[optimal_ncomp]) * 100))
cat(sprintf("  Permutation p-value : %.4f (n = %d)\n", perm_results$p_value, perm_results$n_perm))
cat(sprintf("  VIP > 1.5           : %d metabolites (%.1f%%)\n",
            nrow(important_metabolites),
            nrow(important_metabolites)/nrow(vip_df)*100))
cat(sprintf("  VIP > 2.0           : %d metabolites\n", sum(vip_df$VIP_Score > 2.0)))
cat(sprintf("  VIP > 2.5           : %d metabolites\n", sum(vip_df$VIP_Score > 2.5)))

if (perm_results$p_value < 0.05) {
  cat("✅ Statistical conclusion: Model significantly better than random (p < 0.05)\n")
} else if (perm_results$p_value < 0.1) {
  cat("⚠️  Statistical conclusion: Marginal significance (0.05 < p < 0.10)\n")
} else {
  cat("ℹ️  Statistical conclusion: Not statistically significant (p >= 0.10)\n")
}

cat("\n📁 All results saved in 'output/' directory:\n")
cat("   - VIP_all.csv, VIP_important.csv, VIP_top50.csv\n")
cat("   - permutation_accuracies.csv, permutation_summary.csv\n")
cat("   - PLSDA_performance_metrics.csv\n")
if (optimal_ncomp == 1) cat("   - PLSDA_component1_scores.csv\n")
cat("   - figures/ (VIP, permutation, component selection, ...)\n")

cat("\n🔜 Next step: run additional analyses (e.g., random forest) or generate final report\n")
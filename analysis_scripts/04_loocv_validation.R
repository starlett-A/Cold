# =================================================================
# 🔍 LOOCV VALIDATION (FINAL MODEL ASSESSMENT)
# Purpose: Evaluate PLS-DA performance on important metabolites only
# Input:   PLS-DA results (VIP, optimal ncomp) + scaled data
# Output:  Confusion matrix, accuracy, 95% CI, prediction plots
# =================================================================

cat("🚀 Starting LOOCV validation...\n")
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
# 1. LOAD NECESSARY INPUTS
# =============================================
cat("\n📂 STEP 1: Loading previous results...\n")

# --- Scaled data ---
scaled_file <- "output/preprocessed_data_scaled.csv"
if (!file.exists(scaled_file)) {
  stop("Error: Scaled data not found at: ", scaled_file,
       "\nPlease run 03_pca_analysis.R first.")
}
scaled_data <- read_csv(scaled_file, show_col_types = FALSE)

# --- VIP table (for important metabolite names) ---
vip_file <- "output/VIP_all.csv"
if (!file.exists(vip_file)) {
  stop("Error: VIP table not found at: ", vip_file,
       "\nPlease run 04_plsda_analysis.R first.")
}
vip_df <- read_csv(vip_file, show_col_types = FALSE)

# --- Performance metrics (to get optimal_ncomp) ---
perf_file <- "output/PLSDA_performance_metrics.csv"
if (!file.exists(perf_file)) {
  stop("Error: PLSDA performance metrics not found at: ", perf_file,
       "\nPlease run 04_plsda_analysis.R first.")
}
perf_metrics <- read_csv(perf_file, show_col_types = FALSE)

# Extract optimal ncomp from the saved metrics
optimal_ncomp <- as.numeric(perf_metrics$Value[perf_metrics$Metric == "Optimal components"])
cat(sprintf("  Optimal ncomp = %d (from PLS-DA)\n", optimal_ncomp))

# =============================================
# 2. DEFINE IMPORTANT METABOLITES (VIP > 1.5)
# =============================================
cat("\n🔧 STEP 2: Selecting important metabolites (VIP > 1.5)...\n")

vip_threshold <- 1.5
important_names <- vip_df$Metabolite[vip_df$VIP_Score > vip_threshold]

if (length(important_names) == 0) {
  stop("Error: No metabolites with VIP > ", vip_threshold, " found.")
}
cat(sprintf("  Selected %d important metabolites\n", length(important_names)))

# Ensure they exist in the scaled data
missing <- setdiff(important_names, colnames(scaled_data))
if (length(missing) > 0) {
  cat(sprintf("  ⚠️  %d metabolites missing from data – removed\n", length(missing)))
  important_names <- intersect(important_names, colnames(scaled_data))
}
cat(sprintf("  Final count: %d metabolites\n", length(important_names)))

# =============================================
# 3. PREPARE REDUCED DATA MATRIX
# =============================================
cat("\n📊 STEP 3: Preparing reduced data matrix...\n")

X <- as.matrix(scaled_data[, important_names, drop = FALSE])
rownames(X) <- scaled_data$Sample_ID
Y <- factor(scaled_data$Group, levels = c("Fit-Good", "Fit-Poor"))

cat(sprintf("  Matrix dimensions: %d samples x %d metabolites\n", nrow(X), ncol(X)))

# =============================================
# 4. LEAVE-ONE-OUT CROSS-VALIDATION
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("STEP 4: Leave-One-Out Cross-Validation\n")
cat(paste0(strrep("=", 60), "\n"))

set.seed(2026)
n <- nrow(X)
predictions <- character(n)
true_labels <- as.character(Y)

pb <- txtProgressBar(min = 0, max = n, style = 3)
for (i in 1:n) {
  train_idx <- setdiff(1:n, i)
  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  X_test  <- X[i, , drop = FALSE]
  
  if (length(unique(Y_train)) < 2) {
    predictions[i] <- as.character(Y_train[1])
  } else {
    model <- mixOmics::plsda(X_train, Y_train, ncomp = optimal_ncomp)
    pred <- predict(model, X_test)$class$max.dist[, optimal_ncomp]
    predictions[i] <- as.character(pred)
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# =============================================
# 5. CALCULATE PERFORMANCE METRICS
# =============================================
cat("\n📈 STEP 5: Calculating performance metrics...\n")

# Confusion matrix
cm <- table(Predicted = factor(predictions, levels = levels(Y)),
            True = Y)
accuracy <- sum(diag(cm)) / n

cat(sprintf("  Overall accuracy: %.2f%% (%d/%d)\n", accuracy * 100, sum(diag(cm)), n))
cat("  Confusion matrix:\n")
print(cm)

# Sensitivity / specificity (binary case)
if (nlevels(Y) == 2) {
  TP <- cm["Fit-Good", "Fit-Good"]
  FN <- cm["Fit-Poor", "Fit-Good"]
  TN <- cm["Fit-Poor", "Fit-Poor"]
  FP <- cm["Fit-Good", "Fit-Poor"]
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  balanced_accuracy <- (sensitivity + specificity) / 2
  
  cat(sprintf("  Sensitivity (Fit-Good): %.2f%%\n", sensitivity * 100))
  cat(sprintf("  Specificity (Fit-Poor): %.2f%%\n", specificity * 100))
  cat(sprintf("  Balanced accuracy: %.2f%%\n", balanced_accuracy * 100))
} else {
  sensitivity <- NA
  specificity <- NA
  balanced_accuracy <- NA
}

# 95% binomial confidence interval
binom_test <- binom.test(sum(diag(cm)), n, p = 0.5, alternative = "two.sided")
ci_lower <- binom_test$conf.int[1] * 100
ci_upper <- binom_test$conf.int[2] * 100
cat(sprintf("  95%% CI for accuracy: [%.1f%%, %.1f%%]\n", ci_lower, ci_upper))
cat(sprintf("  Binomial test p-value: %.4f\n", binom_test$p.value))

# =============================================
# 6. SAVE DETAILED RESULTS
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 6: Saving results tables\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

# Detailed per-sample predictions
detailed_df <- data.frame(
  Sample_ID = rownames(X),
  True_Class = true_labels,
  Predicted_Class = predictions,
  Correct = true_labels == predictions
)
write_csv(detailed_df, "output/LOOCV_detailed_predictions.csv")
cat("✅ Saved: output/LOOCV_detailed_predictions.csv\n")

# Performance summary
report_df <- data.frame(
  Metric = c("Optimal_ncomp", "Number_of_metabolites", "Sample_size",
             "Fit_Good_samples", "Fit_Poor_samples", "Overall_accuracy",
             "CI_lower", "CI_upper", "Sensitivity", "Specificity",
             "Balanced_accuracy", "Misclassifications", "Binomial_p_value"),
  Value = c(optimal_ncomp, length(important_names), n,
            sum(Y == "Fit-Good"), sum(Y == "Fit-Poor"),
            accuracy * 100,
            ci_lower, ci_upper,
            ifelse(is.na(sensitivity), NA, sensitivity * 100),
            ifelse(is.na(specificity), NA, specificity * 100),
            ifelse(is.na(balanced_accuracy), NA, balanced_accuracy * 100),
            n - sum(diag(cm)),
            binom_test$p.value)
)
write_csv(report_df, "output/LOOCV_final_report.csv")
cat("✅ Saved: output/LOOCV_final_report.csv\n")

# =============================================
# 7. GENERATE DIAGNOSTIC PLOTS
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 7: Generating LOOCV plots\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# 7.1 Confusion matrix heatmap
cm_df <- as.data.frame(as.table(cm))
colnames(cm_df) <- c("Predicted", "True", "Count")

p_cm <- ggplot(cm_df, aes(x = True, y = Predicted, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Count/n*100)),
            size = 5, color = "black") +
  scale_fill_gradient(low = "white", high = "#2E86AB") +
  labs(title = sprintf("Confusion Matrix (LOOCV)\nAccuracy = %.1f%%", accuracy * 100),
       x = "True Class", y = "Predicted Class") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("output/figures/LOOCV_confusion_matrix.png", p_cm, width = 7, height = 6, dpi = 600)
ggsave("output/figures/LOOCV_confusion_matrix.pdf", p_cm, width = 7, height = 6)

# 7.2 Per-sample prediction accuracy
plot_data <- data.frame(
  Sample = 1:n,
  True_Label = true_labels,
  Predicted_Label = predictions,
  Correct = true_labels == predictions
)

p_pred <- ggplot(plot_data, aes(x = Sample, y = True_Label,
                                color = Correct, shape = Predicted_Label)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "#E15759", "TRUE" = "#4E79A7")) +
  scale_shape_manual(values = c("Fit-Good" = 16, "Fit-Poor" = 17)) +
  labs(title = "LOOCV Predictions per Sample",
       subtitle = sprintf("Accuracy = %.1f%% | n = %d", accuracy * 100, n),
       x = "Sample Index", y = "True Class",
       color = "Correct?", shape = "Predicted") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("output/figures/LOOCV_predictions.png", p_pred, width = 10, height = 5, dpi = 600)
ggsave("output/figures/LOOCV_predictions.pdf", p_pred, width = 10, height = 5)
cat("✅ LOOCV plots saved in output/figures/\n")

# =============================================
# 8. FINAL SUMMARY
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 LOOCV VALIDATION COMPLETE\n")
cat(paste0(strrep("=", 60), "\n\n"))

cat("Model configuration:\n")
cat(sprintf("  Latent components : %d\n", optimal_ncomp))
cat(sprintf("  Important metabolites : %d\n", length(important_names)))
cat(sprintf("  Samples (Fit-Good / Fit-Poor) : %d / %d\n",
            sum(Y == "Fit-Good"), sum(Y == "Fit-Poor")))

cat("\nPerformance:\n")
cat(sprintf("  LOOCV accuracy : %.2f%%\n", accuracy * 100))
cat(sprintf("  95%% CI         : [%.1f%%, %.1f%%]\n", ci_lower, ci_upper))
cat(sprintf("  Binomial p      : %.4f\n", binom_test$p.value))
if (!is.na(balanced_accuracy)) {
  cat(sprintf("  Sensitivity     : %.2f%%\n", sensitivity * 100))
  cat(sprintf("  Specificity     : %.2f%%\n", specificity * 100))
  cat(sprintf("  Balanced accuracy : %.2f%%\n", balanced_accuracy * 100))
}
cat(sprintf("  Misclassifications : %d\n", n - sum(diag(cm))))

if (binom_test$p.value < 0.05) {
  cat("\n✅ Statistical significance: p < 0.05\n")
} else {
  cat("\nℹ️  Statistical significance: p ≥ 0.05 (but may still be meaningful)\n")
}

cat("\n📁 All results saved in 'output/' directory\n")
cat("   - LOOCV_detailed_predictions.csv\n")
cat("   - LOOCV_final_report.csv\n")
cat("   - figures/LOOCV_confusion_matrix.png/pdf\n")
cat("   - figures/LOOCV_predictions.png/pdf\n")
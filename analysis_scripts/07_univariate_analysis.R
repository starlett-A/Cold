# =================================================================
# 📊 UNIVARIATE ANALYSIS (DIFFERENTIAL ABUNDANCE)
# Purpose: Statistical testing + effect size on RF-selected metabolites
# Input:   RF-selected abundance data (log2-transformed)
# Output:  Full results, significant metabolites, summary statistics
# =================================================================

cat("🚀 Starting univariate analysis...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("dplyr", "effsize", "pROC", "ggplot2", "nortest", "MASS",
                "readr", "magrittr", "tidyr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD DATA
# =============================================
cat("\n📂 STEP 1: Loading RF-selected metabolite abundance...\n")

input_file <- "output/RF_selected_metabolites_abundance.csv"
if (!file.exists(input_file)) {
  stop("Error: Abundance file not found at: ", input_file,
       "\nPlease run 07_random_forest.R first.")
}

abundance <- read_csv(input_file, show_col_types = FALSE)

cat(sprintf("  Data loaded: %d samples, %d columns\n",
            nrow(abundance), ncol(abundance)))

# Ensure Group column is factor
if (!"Group" %in% colnames(abundance)) {
  stop("Data must contain a 'Group' column")
}
abundance$Group <- factor(abundance$Group, levels = c("Fit-Good", "Fit-Poor"))

# Identify metabolite columns (everything except Sample_ID and Group)
metabolite_cols <- setdiff(colnames(abundance), c("Sample_ID", "Group"))
cat(sprintf("  Metabolites to test: %d\n", length(metabolite_cols)))

# =============================================
# 2. UNIVARIATE ANALYSIS FUNCTION
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("STEP 2: Running univariate tests (t-test / Wilcoxon)\n")
cat(paste0(strrep("=", 60), "\n"))

univariate_analysis <- function(data, metabolites, group_var) {
  results_list <- list()
  
  pb <- txtProgressBar(min = 0, max = length(metabolites), style = 3)
  
  for (i in seq_along(metabolites)) {
    met <- metabolites[i]
    
    good_vals <- data[[met]][data[[group_var]] == "Fit-Good"]
    poor_vals <- data[[met]][data[[group_var]] == "Fit-Poor"]
    good_vals <- good_vals[!is.na(good_vals)]
    poor_vals <- poor_vals[!is.na(poor_vals)]
    
    if (length(good_vals) < 2 || length(poor_vals) < 2) next
    
    # Normality test
    shapiro_good <- tryCatch(shapiro.test(good_vals)$p.value, error = function(e) NA)
    shapiro_poor <- tryCatch(shapiro.test(poor_vals)$p.value, error = function(e) NA)
    is_normal <- (!is.na(shapiro_good) && !is.na(shapiro_poor) &&
                    shapiro_good > 0.05 && shapiro_poor > 0.05)
    
    # Statistical test
    if (is_normal && length(good_vals) >= 3 && length(poor_vals) >= 3) {
      test_result <- tryCatch(
        t.test(poor_vals, good_vals, var.equal = TRUE),
        error = function(e) wilcox.test(poor_vals, good_vals, exact = FALSE)
      )
      test_type <- "t-test"
    } else {
      test_result <- wilcox.test(poor_vals, good_vals, exact = FALSE)
      test_type <- "Wilcoxon"
    }
    
    # Effect size (Cliff's delta)
    cliff_est <- tryCatch(
      cliff.delta(poor_vals, good_vals)$estimate,
      error = function(e) NA
    )
    
    # Cohen's d (only if normal)
    cohens_d <- NA
    if (is_normal && test_type == "t-test") {
      n1 <- length(good_vals); n2 <- length(poor_vals)
      var1 <- var(good_vals); var2 <- var(poor_vals)
      if (n1 > 1 && n2 > 1 && !is.na(var1) && !is.na(var2) && var1 > 0 && var2 > 0) {
        pooled_sd <- sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1 + n2 - 2))
        if (pooled_sd > 0) cohens_d <- (mean(poor_vals) - mean(good_vals)) / pooled_sd
      }
    }
    
    # Fold change (log2 space → linear)
    log2fc <- mean(poor_vals, na.rm = TRUE) - mean(good_vals, na.rm = TRUE)
    linear_fc <- 2^log2fc
    if (log2fc > 0) {
      direction <- "Upregulated in Poor"
      display_fc <- linear_fc
    } else {
      direction <- "Downregulated in Poor"
      display_fc <- 1 / linear_fc
    }
    
    # AUC
    auc_val <- tryCatch({
      roc_obj <- pROC::roc(response = data[[group_var]],
                           predictor = data[[met]],
                           levels = c("Fit-Good", "Fit-Poor"),
                           direction = "<", quiet = TRUE)
      as.numeric(roc_obj$auc)
    }, error = function(e) NA)
    
    results_list[[met]] <- list(
      Metabolite       = met,
      Sample_Size_Good = length(good_vals),
      Sample_Size_Poor = length(poor_vals),
      Shapiro_Good_p   = shapiro_good,
      Shapiro_Poor_p   = shapiro_poor,
      Normal_Distribution = is_normal,
      Test_Type        = test_type,
      p_value          = test_result$p.value,
      Cliff_Delta      = cliff_est,
      Cohen_d          = cohens_d,
      log2FC           = log2fc,
      Linear_FC        = linear_fc,
      Display_Fold_Change = display_fc,
      Direction        = direction,
      AUC              = auc_val,
      Mean_Good        = mean(good_vals, na.rm = TRUE),
      Mean_Poor        = mean(poor_vals, na.rm = TRUE),
      Median_Good      = median(good_vals, na.rm = TRUE),
      Median_Poor      = median(poor_vals, na.rm = TRUE),
      SD_Good          = sd(good_vals, na.rm = TRUE),
      SD_Poor          = sd(poor_vals, na.rm = TRUE)
    )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Combine and adjust p-values
  results_df <- do.call(rbind, lapply(results_list, as.data.frame, stringsAsFactors = FALSE))
  results_df$p_adj_BH <- p.adjust(results_df$p_value, method = "BH")
  results_df$p_adj_Bonferroni <- p.adjust(results_df$p_value, method = "bonferroni")
  results_df <- results_df[order(results_df$p_value), ]
  results_df$Rank <- seq_len(nrow(results_df))
  
  return(results_df)
}

# Run analysis
univariate_results <- univariate_analysis(abundance, metabolite_cols, "Group")
cat(sprintf("\n✅ Completed testing %d metabolites\n", nrow(univariate_results)))

# =============================================
# 3. DEFINE SIGNIFICANCE CRITERIA & FILTER
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 3: Selecting significant metabolites\n")
cat(paste0(strrep("=", 50), "\n"))

sig_bh    <- univariate_results$p_adj_BH < 0.05 & !is.na(univariate_results$p_adj_BH)
sig_effect <- abs(univariate_results$Cliff_Delta) > 0.5 & !is.na(univariate_results$Cliff_Delta)
sig_fc     <- abs(univariate_results$log2FC) > 0.58 & !is.na(univariate_results$log2FC)

significant_idx <- which(sig_bh | (sig_effect & sig_fc))
n_significant <- length(significant_idx)

cat(sprintf("  FDR < 0.05: %d metabolites\n", sum(sig_bh)))
cat(sprintf("  |Cliff's Delta| > 0.5: %d metabolites\n", sum(sig_effect)))
cat(sprintf("  |log2FC| > 0.58: %d metabolites\n", sum(sig_fc)))
cat(sprintf("  Combined significant: %d metabolites\n", n_significant))

# =============================================
# 4. SAVE RESULTS
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 4: Saving output files\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output")) dir.create("output", recursive = TRUE)

# Full results
write_csv(univariate_results, "output/univariate_all_results.csv")
cat("✅ Saved: output/univariate_all_results.csv\n")

# Significant subset
if (n_significant > 0) {
  sig_df <- univariate_results[significant_idx, ] %>%
    arrange(p_value)
  
  # Key columns
  key_cols <- c("Metabolite", "Rank", "p_value", "p_adj_BH", "p_adj_Bonferroni",
                "log2FC", "Linear_FC", "Display_Fold_Change", "Direction",
                "Cliff_Delta", "Cohen_d", "AUC", "Test_Type",
                "Mean_Good", "Mean_Poor", "Sample_Size_Good", "Sample_Size_Poor")
  available_cols <- intersect(key_cols, names(sig_df))
  sig_df_out <- sig_df[, available_cols]
  
  write_csv(sig_df_out, "output/univariate_significant_metabolites.csv")
  cat(sprintf("✅ Saved: output/univariate_significant_metabolites.csv (%d metabolites)\n", nrow(sig_df_out)))
} else {
  write_csv(data.frame(), "output/univariate_significant_metabolites.csv")
  cat("⚠️  No significant metabolites (empty file written)\n")
}

# Summary table
summary_df <- data.frame(
  Metric = c("Total metabolites",
             "p < 0.05 (raw)",
             "FDR < 0.05",
             "Bonferroni < 0.05",
             "|log2FC| > 0.58",
             "|Cliff's Delta| > 0.5",
             "Significant (combined)"),
  Count = c(nrow(univariate_results),
            sum(univariate_results$p_value < 0.05, na.rm = TRUE),
            sum(sig_bh),
            sum(univariate_results$p_adj_Bonferroni < 0.05, na.rm = TRUE),
            sum(sig_fc),
            sum(sig_effect),
            n_significant)
)
write_csv(summary_df, "output/univariate_summary.csv")
cat("✅ Saved: output/univariate_summary.csv\n")

# =============================================
# 5. FINAL REPORT
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 UNIVARIATE ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 60), "\n\n"))

cat("Statistical tests used:\n")
print(table(univariate_results$Test_Type))

cat(sprintf("\n  Total tested: %d\n", nrow(univariate_results)))
cat(sprintf("  Significant : %d\n", n_significant))
if (n_significant > 0) {
  n_up <- sum(univariate_results$log2FC[significant_idx] > 0)
  n_down <- n_significant - n_up
  cat(sprintf("  Up in Poor   : %d\n", n_up))
  cat(sprintf("  Down in Poor : %d\n", n_down))
  
  cat("\n🏆 Top 5 most significant metabolites:\n")
  top5 <- head(sig_df, 5)
  for (i in 1:nrow(top5)) {
    cat(sprintf("  %d. %-30s log2FC=%.2f  p=%.2e  FDR=%.2e  Cliff=%.2f\n",
                i, top5$Metabolite[i], top5$log2FC[i], top5$p_value[i],
                top5$p_adj_BH[i], top5$Cliff_Delta[i]))
  }
}

cat("\n📁 Key outputs:\n")
cat("  output/univariate_all_results.csv\n")
cat("  output/univariate_significant_metabolites.csv\n")
cat("  output/univariate_summary.csv\n")

cat("\n🔜 Next step: functional enrichment, integration with 16S, or figure generation.\n")
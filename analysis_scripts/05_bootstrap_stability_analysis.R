# =================================================================
# 🔄 BOOTSTRAP STABILITY ANALYSIS
# Purpose: Evaluate selection stability of important metabolites
# Input:   Scaled data, VIP list, optimal ncomp, name mapping
# Output:  Frequency tables, confidence levels, diagnostic plots
# =================================================================

cat("🚀 Starting Bootstrap stability analysis...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("dplyr", "ggplot2", "ggrepel", "mixOmics", "tidyr", "readr", "magrittr"))
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

# VIP table (important metabolites)
vip_file <- "output/VIP_important.csv"
if (!file.exists(vip_file)) {
  stop("Error: VIP important table not found at: ", vip_file,
       "\nPlease run 04_plsda_analysis.R first.")
}
vip_df <- read_csv(vip_file, show_col_types = FALSE)

# Optimal ncomp (from PLSDA performance metrics)
perf_file <- "output/PLSDA_performance_metrics.csv"
if (!file.exists(perf_file)) {
  stop("Error: PLSDA performance metrics not found at: ", perf_file,
       "\nPlease run 04_plsda_analysis.R first.")
}
perf_metrics <- read_csv(perf_file, show_col_types = FALSE)
optimal_ncomp <- as.numeric(perf_metrics$Value[perf_metrics$Metric == "Optimal components"])
cat(sprintf("  Optimal ncomp = %d\n", optimal_ncomp))

# Name mapping (for original metabolite labels)
mapping_file <- "output/metabolite_name_mapping.csv"
name_mapping <- NULL
if (file.exists(mapping_file)) {
  name_mapping <- read_csv(mapping_file, show_col_types = FALSE)
  cat("✅ Name mapping loaded\n")
} else {
  cat("⚠️  Name mapping not found, using safe names\n")
}

# =============================================
# 2. PREPARE DATA & METABOLITE LIST
# =============================================
cat("\n🔧 STEP 2: Preparing data matrix and metabolite list...\n")

# Data matrix
X <- as.matrix(scaled_data[, 3:ncol(scaled_data)])
rownames(X) <- scaled_data$Sample_ID
Y <- factor(scaled_data$Group, levels = c("Fit-Good", "Fit-Poor"))

# Important metabolite names (safe names)
metabolite_list <- vip_df$Metabolite
if (is.null(metabolite_list)) stop("VIP table does not contain 'Metabolite' column")
cat(sprintf("  Candidate metabolites from VIP: %d\n", length(metabolite_list)))

# Ensure they are in the data
missing <- setdiff(metabolite_list, colnames(X))
if (length(missing) > 0) {
  cat(sprintf("  ⚠️  %d metabolites missing from data – removed\n", length(missing)))
  metabolite_list <- intersect(metabolite_list, colnames(X))
}
cat(sprintf("  Final list: %d metabolites\n", length(metabolite_list)))

# =============================================
# 3. BOOTSTRAP STABILITY FUNCTION
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("STEP 3: Bootstrap sampling (n = 100)\n")
cat(paste0(strrep("=", 60), "\n"))

bootstrap_stability <- function(X, Y, ncomp, metabolite_list,
                                n_bootstrap = 100, vip_threshold = 1.5,
                                seed = 2026) {
  set.seed(seed)
  n_samples <- nrow(X)
  
  # Initialize counters
  selection_count <- setNames(rep(0, length(metabolite_list)), metabolite_list)
  successful <- 0
  failed <- 0
  
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  for (b in 1:n_bootstrap) {
    # Resample with replacement
    idx <- sample(1:n_samples, size = n_samples, replace = TRUE)
    X_boot <- X[idx, , drop = FALSE]
    Y_boot <- Y[idx]
    
    # Skip if only one class
    if (length(unique(Y_boot)) < 2) {
      failed <- failed + 1
      setTxtProgressBar(pb, b)
      next
    }
    
    # PLS-DA on bootstrap sample
    model_boot <- tryCatch({
      mixOmics::plsda(X_boot, Y_boot, ncomp = ncomp)
    }, error = function(e) NULL)
    
    if (is.null(model_boot)) {
      failed <- failed + 1
      setTxtProgressBar(pb, b)
      next
    }
    
    # VIP scores
    vip_boot <- tryCatch({
      vip_mat <- mixOmics::vip(model_boot)
      if (ncomp > 1) apply(vip_mat[, 1:ncomp, drop = FALSE], 1, mean)
      else as.numeric(vip_mat[, 1])
    }, error = function(e) NULL)
    
    if (is.null(vip_boot)) {
      failed <- failed + 1
      setTxtProgressBar(pb, b)
      next
    }
    names(vip_boot) <- rownames(vip(model_boot))
    
    # Selected metabolites in this iteration
    selected <- names(vip_boot)[vip_boot > vip_threshold]
    found <- intersect(metabolite_list, selected)
    if (length(found) > 0) selection_count[found] <- selection_count[found] + 1
    successful <- successful + 1
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  if (successful == 0) stop("All bootstrap iterations failed. Check data or parameters.")
  
  freq <- (selection_count / successful) * 100
  result_df <- data.frame(
    Metabolite = names(freq),
    Bootstrap_Frequency = round(freq, 2),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Bootstrap_Frequency)) %>%
    mutate(
      Rank = row_number(),
      Frequency_Category = case_when(
        Bootstrap_Frequency >= 75 ~ "High (≥75%)",
        Bootstrap_Frequency >= 50 ~ "Medium (50-75%)",
        Bootstrap_Frequency >= 25 ~ "Low (25-50%)",
        TRUE ~ "Very Low (<25%)"
      )
    )
  
  list(
    frequency_table = result_df,
    selection_count = selection_count,
    successful = successful,
    total = n_bootstrap,
    failed = failed
  )
}

# Run bootstrap
boot_start <- Sys.time()
boot_results <- bootstrap_stability(
  X = X, Y = Y,
  ncomp = optimal_ncomp,
  metabolite_list = metabolite_list,
  n_bootstrap = 100,
  vip_threshold = 1.5,
  seed = 2026
)
boot_end <- Sys.time()
cat(sprintf("\n⏱️  Bootstrap runtime: %.1f min\n",
            difftime(boot_end, boot_start, units = "mins")))

# =============================================
# 4. MAP NAMES AND CREATE FINAL TABLE
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 4: Integrating VIP & Bootstrap results\n")
cat(paste0(strrep("=", 50), "\n"))

# Scale function for composite score
scale01 <- function(x) {
  if (length(unique(na.omit(x))) <= 1) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Merge with VIP original scores
boot_df <- boot_results$frequency_table
vip_sub <- vip_df %>% select(Metabolite, VIP_Score)
final_tab <- boot_df %>%
  left_join(vip_sub, by = "Metabolite")

# Add original names if mapping exists
if (!is.null(name_mapping)) {
  final_tab <- final_tab %>%
    left_join(name_mapping %>% select(Safe_Name, Original_Name),
              by = c("Metabolite" = "Safe_Name")) %>%
    mutate(Display_Name = ifelse(is.na(Original_Name), Metabolite, Original_Name))
} else {
  final_tab <- final_tab %>% mutate(Display_Name = Metabolite)
}

# VIP normalized and composite score
final_tab <- final_tab %>%
  mutate(
    VIP_Normalized = scale01(VIP_Score),
    Freq_Normalized = Bootstrap_Frequency / 100,
    Composite_Score = 0.6 * VIP_Normalized + 0.4 * Freq_Normalized,
    Confidence_Level = case_when(
      Bootstrap_Frequency >= 75 & VIP_Score >= 2.0 ~ "Very High",
      Bootstrap_Frequency >= 50 & VIP_Score >= 1.5 ~ "High",
      Bootstrap_Frequency >= 25 | VIP_Score >= 1.5 ~ "Medium",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(desc(Composite_Score)) %>%
  mutate(Rank = row_number())

# Save main table
write_csv(final_tab, "output/bootstrap_stability_results.csv")
cat("✅ Saved: output/bootstrap_stability_results.csv\n")

# =============================================
# 5. SAVE HIGH-CONFIDENCE SUBSETS
# =============================================
cat("\n💾 Saving high-confidence subsets...\n")
for (level in c("Very High", "High")) {
  sub <- final_tab %>% filter(Confidence_Level == level)
  if (nrow(sub) > 0) {
    fname <- paste0("output/bootstrap_", gsub(" ", "_", tolower(level)), "_confidence.csv")
    write_csv(sub, fname)
    cat(sprintf("  %s confidence: %d metabolites → %s\n", level, nrow(sub), fname))
  }
}

# =============================================
# 6. VISUALIZATIONS
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 5: Generating diagnostic plots\n")
cat(paste0(strrep("=", 50), "\n"))

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# 6.1 Frequency distribution histogram
p_freq <- ggplot(final_tab, aes(x = Bootstrap_Frequency)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = c(25, 50, 75), linetype = "dashed",
             color = c("orange", "red", "darkred")) +
  annotate("text", x = c(25, 50, 75), y = Inf,
           label = c("25%", "50%", "75%"), vjust = 1.5, hjust = -0.2,
           color = c("orange", "red", "darkred")) +
  labs(title = "Bootstrap Frequency Distribution",
       subtitle = sprintf("Mean frequency: %.1f%%", mean(final_tab$Bootstrap_Frequency)),
       x = "Bootstrap Frequency (%)", y = "Number of metabolites") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

ggsave("output/figures/bootstrap_frequency_distribution.png", p_freq, width = 10, height = 6, dpi = 600)
ggsave("output/figures/bootstrap_frequency_distribution.pdf", p_freq, width = 10, height = 6)

# 6.2 Top 20 barplot
top_n <- min(20, nrow(final_tab))
top_metabs <- final_tab[1:top_n, ]
top_metabs$Display_Name <- paste0(top_metabs$Rank, ". ", substr(top_metabs$Display_Name, 1, 30))
top_metabs$Display_Name <- factor(top_metabs$Display_Name,
                                  levels = top_metabs$Display_Name[order(top_metabs$Bootstrap_Frequency)])

p_bar <- ggplot(top_metabs, aes(x = Bootstrap_Frequency, y = Display_Name, fill = Confidence_Level)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Bootstrap_Frequency)), hjust = -0.1, size = 3.5) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Very High" = "darkred", "High" = "red",
                               "Medium" = "orange", "Low" = "gray")) +
  labs(title = paste0("Top ", top_n, " Metabolites by Bootstrap Stability"),
       x = "Bootstrap Frequency (%)", y = NULL,
       fill = "Confidence") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

ggsave("output/figures/bootstrap_top20_barplot.png", p_bar, width = 14, height = 10, dpi = 600)
ggsave("output/figures/bootstrap_top20_barplot.pdf", p_bar, width = 14, height = 10)

# 6.3 VIP vs Bootstrap scatter (only if VIP is available)
if (all(!is.na(final_tab$VIP_Score))) {
  cor_test <- cor.test(final_tab$VIP_Score, final_tab$Bootstrap_Frequency)
  top_label_n <- min(15, nrow(final_tab))
  
  p_scatter <- ggplot(final_tab, aes(x = VIP_Score, y = Bootstrap_Frequency)) +
    geom_point(aes(color = Confidence_Level, size = Composite_Score), alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
    geom_text_repel(data = final_tab[1:top_label_n, ],
                    aes(label = paste0(Rank, ". ", substr(Display_Name, 1, 25))),
                    size = 3, box.padding = 0.5, max.overlaps = 20) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "blue") +
    labs(title = "VIP Score vs Bootstrap Stability",
         subtitle = sprintf("Pearson r = %.2f (p = %.3f)", cor_test$estimate, cor_test$p.value),
         x = "VIP Score", y = "Bootstrap Frequency (%)",
         color = "Confidence", size = "Composite") +
    scale_color_manual(values = c("Very High" = "darkred", "High" = "red",
                                  "Medium" = "orange", "Low" = "gray")) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave("output/figures/bootstrap_VIP_vs_frequency.png", p_scatter, width = 12, height = 8, dpi = 600)
  ggsave("output/figures/bootstrap_VIP_vs_frequency.pdf", p_scatter, width = 12, height = 8)
}
cat("✅ All diagnostic plots saved\n")

# =============================================
# 7. FINAL SUMMARY
# =============================================
cat(paste0("\n", strrep("=", 60), "\n"))
cat("🎉 BOOTSTRAP STABILITY ANALYSIS COMPLETE\n")
cat(paste0(strrep("=", 60), "\n\n"))

cat("Bootstrap settings:\n")
cat(sprintf("  Iterations: %d (successful: %d, failed: %d)\n",
            boot_results$total, boot_results$successful, boot_results$failed))
cat(sprintf("  VIP threshold: > 1.5\n"))

cat("\nStability statistics:\n")
cat(sprintf("  Mean frequency: %.1f%%\n", mean(final_tab$Bootstrap_Frequency)))
cat(sprintf("  ≥ 75%% : %d metabolites\n", sum(final_tab$Bootstrap_Frequency >= 75)))
cat(sprintf("  ≥ 50%% : %d metabolites\n", sum(final_tab$Bootstrap_Frequency >= 50)))
cat(sprintf("  ≥ 25%% : %d metabolites\n", sum(final_tab$Bootstrap_Frequency >= 25)))

cat("\nTop 5 most stable metabolites:\n")
for (i in 1:min(5, nrow(final_tab))) {
  cat(sprintf("  %d. %s (VIP=%.2f, Freq=%.1f%%)\n",
              i, final_tab$Display_Name[i], final_tab$VIP_Score[i],
              final_tab$Bootstrap_Frequency[i]))
}

cat("\n📁 Key output files:\n")
cat("  output/bootstrap_stability_results.csv\n")
cat("  output/bootstrap_high_confidence.csv\n")
cat("  output/figures/ (distribution, top20, scatter)\n")

cat("\n🔜 Next steps: pathway enrichment of highly stable metabolites, biological validation.\n")
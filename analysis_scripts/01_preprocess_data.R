# =================================================================
# 🔄 METABOLOMICS DATA PREPROCESSING PIPELINE
# Purpose: Complete preprocessing workflow for metabolomics data
# Input: Raw metabolomics data (CSV format)
# Output: Preprocessing data ready for statistical analysis
# =================================================================

cat("🚀 Starting metabolomics data preprocessing...\n")
cat("=============================================")

# =============================================
# 1. LOAD AND INSPECT DATA
# =============================================

cat("\n📂 STEP 1: Loading data...\n")

# Use relative path for portability
input_file <- "data/raw_metabolomics_data.csv"

# Check if file exists
if (!file.exists(input_file)){
  stop("Error: Input file not found at:",input_file,
       "\nPlease ensure your data file is in the 'data/' folder.")
}

# 1.1 Load raw data
combined_data <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# 1.2 Basic data information
cat("Data overview:\n")
cat(springf(" Total samples: %d\n", nrow(combined_data)))
cat(springf(" Total variables: %d\n", ncol(combined_data)))

# 1.3 Check group information
cat("\nGroup distribution:\n")
group_counts <- table(combined_data$Group)
print(group_counts)
cat(sprintf("  Fit-Good group: %d samples\n", group_counts["Fit-Good"]))
cat(sprintf("  Fit-Poor group: %d samples\n", group_counts["Fit-Poor"]))

# 1.4 Count metabolites
n_metabolites <- ncol(combined_data) - 2
cat(sprintf("\nNumber of metabolites: %d\n", n_metabolites))

# 1.5 Create valid metabolite names
cat("\nProcessing metabolite names...\n")
original_names <- colnames(combined_data)[3:ncol(combined_data)]
safe_names <- make.names(original_names, unique = TRUE)

# Save name mapping
name_mapping <- data.frame(
  Original_Name = original_names,
  Safe_Name = safe_names,
  stringsAsFactors = FALSE
)

# Save mapping table
if (!dir.exists("output")) dir.create("output", recursive = TRUE)
write.csv(name_mapping, "output/metabolite_name_mapping.csv", row.names = FALSE)
cat("✅ Name mapping saved: output/metabolite_name_mapping.csv\n")

# Apply safe names to data
colnames(combined_data)[1] <- "Sample_ID"
colnames(combined_data)[2] <- "Group"
colnames(combined_data)[3:ncol(combined_data)] <- safe_names

# ============================================
# 2. DATA QUALITY ASSESSMENT
# ============================================

cat("\n🔍 STEP 2: Data quality assessment...\n")

# Extract numeric data
numeric_data <- combined_data[, 3:ncol(combined_data)]
missing_by_sample <- rowSums(is.na(numeric_data))
missing_by_metabolite <- colSums(is.na(numeric_data))

cat("\nMissing value statistics:\n")
overall_missing <- sum(is.na(numeric_data)) / (nrow(numeric_data) * ncol(numeric_data)) * 100
cat(sprintf("  Overall missing rate: %.2f%%\n", overall_missing))
cat(sprintf("  Median missing per sample: %.1f metabolites\n", median(missing_by_sample)))
cat(sprintf("  Median missing per metabolite: %.1f samples\n", median(missing_by_metabolite)))

# ============================================
# 3. MISSING VALUE FILTERING
# ============================================

cat("\n🧹 STEP 3: Missing value filtering...\n")

# 3.1 Set filtering threshold (keep metabolites present in ≥50% of samples)
threshold <- 0.5
cat(sprintf("\nFiltering threshold: Keep metabolites present in ≥%.0f%% of samples\n", threshold * 100))

# 3.2 Calculate presence rate for each metabolite
presence_rate <- apply(numeric_data, 2, function(x) {
  sum(!is.na(x)) / length(x)
})

# 3.3 Filter metabolites
metabolites_to_keep <- names(presence_rate)[presence_rate >= threshold]
metabolites_removed <- names(presence_rate)[presence_rate < threshold]

cat(sprintf("  Original metabolites: %d\n", n_metabolites))
cat(sprintf("  After filtering: %d metabolites\n", length(metabolites_to_keep)))
cat(sprintf("  Removed: %d metabolites\n", length(metabolites_removed)))

# 3.4 Save removed metabolites information
if (length(metabolites_removed) > 0) {
  removed_metabolites_info <- data.frame(
    Metabolite = metabolites_removed,
    Presence_Rate = presence_rate[metabolites_removed],
    Reason = "Missing rate too high",
    stringsAsFactors = FALSE
  )
  write.csv(removed_metabolites_info, "output/removed_metabolites.csv", row.names = FALSE)
  cat("✅ Removed metabolites list saved\n")
} else {
  # Create empty file if no metabolites were removed
  removed_metabolites_info <- data.frame(
    Metabolite = character(),
    Presence_Rate = numeric(),
    Reason = character(),
    stringsAsFactors = FALSE
  )
  write.csv(removed_metabolites_info, "output/removed_metabolites.csv", row.names = FALSE)
  cat("✅ All metabolites passed filtering (empty list created)\n")
}

# 3.5 Create filtered dataset
filtered_data <- data.frame(
  Sample_ID = combined_data$Sample_ID,
  Group = combined_data$Group,
  numeric_data[, metabolites_to_keep, drop = FALSE],
  stringsAsFactors = FALSE
)

# ============================================
# 4. MISSING VALUE IMPUTATION
# ============================================

cat("\n🔄 STEP 4: Missing value imputation...\n")

# 4.1 Check missing values before imputation
numeric_filtered <- filtered_data[, 3:ncol(filtered_data)]
missing_before <- sum(is.na(numeric_filtered)) / length(numeric_filtered) * 100
cat(sprintf("Missing rate before imputation: %.2f%%\n", missing_before))

# 4.2 Impute with half of minimum value (suitable for small datasets)
imputed_data <- numeric_filtered

for (i in 1:ncol(imputed_data)) {
  col_data <- imputed_data[, i]
  non_missing <- col_data[!is.na(col_data)]
  
  if (length(non_missing) > 0) {
    # Find minimum positive value
    positive_vals <- non_missing[non_missing > 0]
    if (length(positive_vals) > 0) {
      min_value <- min(positive_vals, na.rm = TRUE)
      impute_value <- min_value * 0.5
    } else {
      # If no positive values, use overall minimum
      min_value <- min(non_missing, na.rm = TRUE)
      impute_value <- min_value * 0.5
    }
    
    # Impute missing values
    imputed_data[is.na(imputed_data[, i]), i] <- impute_value
  }
}

# 4.3 Verify imputation
missing_after <- sum(is.na(imputed_data)) / length(imputed_data) * 100
cat(sprintf("Missing rate after imputation: %.2f%%\n", missing_after))

# 4.4 Update dataframe
imputed_df <- data.frame(
  Sample_ID = filtered_data$Sample_ID,
  Group = filtered_data$Group,
  imputed_data,
  stringsAsFactors = FALSE
)

# ============================================
# 5. PQN NORMALIZATION
# ============================================

cat("\n⚖️ STEP 5: PQN normalization...\n")

# 5.1 Define PQN normalization function
pqn_normalize <- function(data_matrix) {
  # Ensure no zeros or negative values
  data_matrix[data_matrix <= 0] <- min(data_matrix[data_matrix > 0]) * 0.01
  
  # Calculate reference spectrum (median of all samples)
  reference_spectrum <- apply(data_matrix, 2, median, na.rm = TRUE)
  
  # Calculate quotient matrix
  quotient_matrix <- sweep(data_matrix, 2, reference_spectrum, FUN = "/")
  
  # Calculate correction factors (median of quotients for each sample)
  correction_factors <- apply(quotient_matrix, 1, median, na.rm = TRUE)
  
  # Apply correction factors
  normalized_matrix <- sweep(data_matrix, 1, correction_factors, FUN = "/")
  
  return(normalized_matrix)
}

# 5.2 Prepare data matrix
data_matrix <- as.matrix(imputed_df[, 3:ncol(imputed_df)])
rownames(data_matrix) <- imputed_df$Sample_ID

# 5.3 Perform PQN normalization
cat("Performing PQN normalization...\n")
pqn_matrix <- pqn_normalize(data_matrix)

# 5.4 Check normalization effect
cat("Normalization quality check:\n")
cat(sprintf("  Raw data median: %.4f\n", median(data_matrix)))
cat(sprintf("  Normalized data median: %.4f\n", median(pqn_matrix)))

# 5.5 Create normalized dataframe
normalized_data <- data.frame(
  Sample_ID = imputed_df$Sample_ID,
  Group = imputed_df$Group,
  pqn_matrix,
  stringsAsFactors = FALSE
)

# ============================================
# 6. LOG2 TRANSFORMATION
# ============================================

cat("\n📊 STEP 6: Log2 transformation...\n")

# 6.1 Apply log2(x+1) transformation
log_transformed <- normalized_data
log_transformed[, 3:ncol(log_transformed)] <- log2(normalized_data[, 3:ncol(normalized_data)] + 1)

# 6.2 Check transformation
cat("Data distribution after log2 transformation:\n")
summary_stats <- summary(as.vector(as.matrix(log_transformed[, 3:min(8, ncol(log_transformed))])))
print(summary_stats)

# ============================================
# 7. SAVE PREPROCESSED DATA
# ============================================

cat("\n💾 STEP 7: Saving preprocessed data...\n")

# 7.1 Save preprocessed data
write.csv(log_transformed, 
          "output/preprocessed_data.csv", 
          row.names = FALSE)
cat("✅ Preprocessed data saved: output/preprocessed_data.csv\n")

# 7.2 Save matrix format
matrix_data <- as.matrix(log_transformed[, 3:ncol(log_transformed)])
rownames(matrix_data) <- log_transformed$Sample_ID
write.csv(matrix_data, "output/preprocessed_matrix.csv")
cat("✅ Matrix format saved: output/preprocessed_matrix.csv\n")

# 7.3 Save group information
group_info <- data.frame(
  Sample_ID = log_transformed$Sample_ID,
  Group = log_transformed$Group
)
write.csv(group_info, "output/sample_group_info.csv", row.names = FALSE)
cat("✅ Group information saved: output/sample_group_info.csv\n")

# 7.4 Generate preprocessing report
preprocessing_report <- data.frame(
  Step = c("Raw data", "After filtering", "After imputation", 
           "After PQN normalization", "After log2 transformation"),
  Samples = rep(nrow(combined_data), 5),
  Metabolites = c(
    ncol(combined_data) - 2,
    length(metabolites_to_keep),
    length(metabolites_to_keep),
    length(metabolites_to_keep),
    length(metabolites_to_keep)
  ),
  Description = c(
    "Raw data loaded from CSV",
    paste0("Removed metabolites with missing rate > ", (1-threshold)*100, "%"),
    "Imputed remaining missing values with half-minimum",
    "Applied Probabilistic Quotient Normalization (PQN)",
    "Applied log2(x+1) transformation"
  )
)

write.csv(preprocessing_report, "output/preprocessing_report.csv", row.names = FALSE)
cat("✅ Preprocessing report saved: output/preprocessing_report.csv\n")

# ============================================
# 8. DATA VISUALIZATION
# ============================================

cat("\n📈 STEP 8: Generating quality control plots...\n")

# 8.1 Boxplot of data distribution
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# Create PDF with boxplots
pdf("output/figures/data_distribution_boxplots.pdf", width = 14, height = 8)

par(mfrow = c(1, 2))

# Original data (first 50 metabolites)
boxplot(data_matrix[, 1:min(50, ncol(data_matrix))], 
        main = "Original data distribution\n(Top 50 metabolites)",
        col = "lightblue",
        outline = FALSE,
        las = 2,
        cex.axis = 0.7,
        ylab = "Intensity")

# Normalized data (first 50 metabolites)
boxplot(pqn_matrix[, 1:min(50, ncol(pqn_matrix))], 
        main = "PQN normalized data distribution\n(Top 50 metabolites)",
        col = "lightgreen",
        outline = FALSE,
        las = 2,
        cex.axis = 0.7,
        ylab = "Normalized intensity")

dev.off()
cat("✅ QC plots saved: output/figures/data_distribution_boxplots.pdf\n")

# ============================================
# 9. FINAL SUMMARY
# ============================================

cat("\n" + strrep("=", 50) + "\n")
cat("✅ PREPROCESSING COMPLETE - FINAL SUMMARY\n")
cat(strrep("=", 50) + "\n\n")

cat("📊 FINAL DATA INFORMATION:\n")
cat(sprintf("  Total samples: %d\n", nrow(log_transformed)))
cat(sprintf("  Fit-Good group: %d samples\n", sum(log_transformed$Group == "Fit-Good")))
cat(sprintf("  Fit-Poor group: %d samples\n", sum(log_transformed$Group == "Fit-Poor")))
cat(sprintf("  Final metabolites: %d\n", ncol(log_transformed) - 2))

# Data integrity check
cat("\n🔍 DATA INTEGRITY CHECK:\n")
final_data <- log_transformed[, 3:ncol(log_transformed)]
missing_final <- sum(is.na(final_data)) / length(final_data) * 100
cat(sprintf("  Final missing rate: %.2f%%\n", missing_final))
cat(sprintf("  Minimum value: %.4f\n", min(final_data, na.rm = TRUE)))
cat(sprintf("  Maximum value: %.4f\n", max(final_data, na.rm = TRUE)))
cat(sprintf("  Median value: %.4f\n", median(as.matrix(final_data), na.rm = TRUE)))

cat("\n📁 OUTPUT FILES GENERATED:\n")
cat("  1. output/preprocessed_data.csv - Complete preprocessed dataset\n")
cat("  2. output/preprocessed_matrix.csv - Matrix format for analysis\n")
cat("  3. output/sample_group_info.csv - Sample grouping information\n")
cat("  4. output/metabolite_name_mapping.csv - Metabolite name mapping\n")
cat("  5. output/removed_metabolites.csv - List of removed metabolites\n")
cat("  6. output/preprocessing_report.csv - Step-by-step processing report\n")
cat("  7. output/figures/data_distribution_boxplots.pdf - QC visualization\n")

cat("\n" + strrep("=", 50) + "\n")
cat("🎉 PREPROCESSING PIPELINE FINISHED SUCCESSFULLY!\n")
cat(strrep("=", 50) + "\n")
cat("Data is now ready for statistical analysis.\n")
cat("Next step: Run statistical tests and generate figures.\n")


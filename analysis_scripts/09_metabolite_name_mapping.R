# =================================================================
# 🔄 METABOLITE NAME MAPPING (SAFE → ORIGINAL)
# Purpose: Replace safe R names with original metabolite names in
#          all key output files for biological interpretation.
# Input:   Main output CSVs + metabolite_name_mapping.csv
# Output:  *_Mapped.csv files with original names
# =================================================================

cat("🚀 Starting metabolite name mapping...\n")
cat(paste0(strrep("=", 45), "\n"))

# =============================================
# 0. LOAD REQUIRED PACKAGES
# =============================================
cat("\n📦 Loading required packages...\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c("dplyr", "readr"))
cat("✅ Package loading complete\n")

# =============================================
# 1. LOAD NAME MAPPING FILE
# =============================================
cat("\n📂 STEP 1: Loading metabolite name mapping...\n")

mapping_file <- "output/metabolite_name_mapping.csv"
if (!file.exists(mapping_file)) {
  stop("Error: Name mapping file not found at: ", mapping_file,
       "\nPlease run 02_preprocessing_pipeline.R first to generate it.")
}

name_map <- read_csv(mapping_file, show_col_types = FALSE)
cat(sprintf("  Loaded %d name mappings\n", nrow(name_map)))

# Ensure correct column names
if (!"Safe_Name" %in% colnames(name_map) || !"Original_Name" %in% colnames(name_map)) {
  # Assume first two columns are Safe_Name, Original_Name
  if (ncol(name_map) >= 2) {
    colnames(name_map)[1:2] <- c("Safe_Name", "Original_Name")
    cat("  Renamed first two columns to Safe_Name / Original_Name\n")
  } else {
    stop("Mapping file must contain at least two columns.")
  }
}

# =============================================
# 2. MAPPING HELPER FUNCTION
# =============================================

# Function to replace safe names with original names (using base R for stability)
map_safe_to_original <- function(safe_names) {
  idx <- match(safe_names, name_map$Safe_Name)
  original <- ifelse(is.na(idx), safe_names, name_map$Original_Name[idx])
  return(original)
}

# General mapping function for different file types
process_file <- function(file_path) {
  cat(sprintf("\n📄 Processing: %s\n", basename(file_path)))
  
  if (!file.exists(file_path)) {
    cat("  ⚠️  File not found, skipping.\n")
    return(invisible(NULL))
  }
  
  # Read data (preserving column names as-is)
  data <- read_csv(file_path, show_col_types = FALSE, name_repair = "minimal")
  cat(sprintf("  Dimensions: %d rows x %d columns\n", nrow(data), ncol(data)))
  
  changed <- 0
  
  # CASE 1: File contains a 'Metabolite' column (e.g., VIP tables, RF importance)
  if ("Metabolite" %in% colnames(data)) {
    cat("  Found 'Metabolite' column → mapping values...\n")
    old_vals <- data$Metabolite
    data$Metabolite <- map_safe_to_original(old_vals)
    changed <- sum(old_vals != data$Metabolite, na.rm = TRUE)
    cat(sprintf("  %d / %d metabolite names mapped\n", changed, length(old_vals)))
  }
  
  # CASE 2: Abundance data (columns are metabolite names, first two are Sample_ID, Group)
  if ("Sample_ID" %in% colnames(data) && "Group" %in% colnames(data)) {
    # Identify metabolite columns (those not in the reserved list)
    reserved <- c("Sample_ID", "Group")
    met_cols <- setdiff(colnames(data), reserved)
    
    if (length(met_cols) > 0) {
      # Only map if these look like safe names (e.g., start with X or contain dots)
      # We'll attempt mapping; unchanged columns will keep original name
      new_names <- map_safe_to_original(met_cols)
      idx_changed <- which(met_cols != new_names)
      if (length(idx_changed) > 0) {
        old_names <- met_cols[idx_changed]
        new_names_sub <- new_names[idx_changed]
        # Update column names
        colnames(data)[which(colnames(data) %in% old_names)] <- new_names_sub
        changed <- length(idx_changed)
        cat(sprintf("  %d column names mapped (abundance matrix)\n", changed))
      } else {
        cat("  No column names needed mapping.\n")
      }
    }
  }
  
  # CASE 3: Generic – search for any character column that contains safe names
  if (changed == 0 && !("Metabolite" %in% colnames(data))) {
    char_cols <- sapply(data, is.character)
    for (col in names(data)[char_cols]) {
      if (any(data[[col]] %in% name_map$Safe_Name, na.rm = TRUE)) {
        cat(sprintf("  Found character column '%s' with safe names → mapping...\n", col))
        old_vals <- data[[col]]
        data[[col]] <- map_safe_to_original(old_vals)
        changed <- sum(old_vals != data[[col]], na.rm = TRUE)
        cat(sprintf("  %d values mapped in column '%s'\n", changed, col))
        break  # only map the first matched column
      }
    }
  }
  
  # Save mapped version
  out_path <- gsub("\\.csv$", "_Mapped.csv", file_path)
  write_csv(data, out_path)
  cat(sprintf("  ✅ Saved: %s\n", basename(out_path)))
  
  return(invisible(data))
}

# =============================================
# 3. APPLY TO KEY OUTPUT FILES
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 2: Mapping metabolite names in output files\n")
cat(paste0(strrep("=", 50), "\n"))

# List of output files that contain metabolite names
files_to_map <- c(
  "output/RF_selected_metabolites_abundance.csv",
  "output/univariate_all_results.csv",
  "output/univariate_significant_metabolites.csv",
  "output/VIP_all.csv",
  "output/VIP_important.csv",
  "output/VIP_top50.csv",
  "output/bootstrap_stability_results.csv",
  "output/RandomForest_importance.csv"
)

# Process files that exist
for (f in files_to_map) {
  process_file(f)
}

# =============================================
# 4. SUMMARY OF MAPPING COVERAGE (OPTIONAL)
# =============================================
cat(paste0("\n", strrep("=", 50), "\n"))
cat("STEP 3: Mapping coverage summary\n")
cat(paste0(strrep("=", 50), "\n"))

# Use the univariate results (if exists) to show how many unique metabolites could be mapped
uni_file <- "output/univariate_all_results.csv"
if (file.exists(uni_file)) {
  uni_data <- read_csv(uni_file, show_col_types = FALSE)
  if ("Metabolite" %in% colnames(uni_data)) {
    all_safe <- unique(uni_data$Metabolite)
    mappable <- intersect(all_safe, name_map$Safe_Name)
    cat(sprintf("  Total unique metabolites in univariate results: %d\n", length(all_safe)))
    cat(sprintf("  Successfully mapped: %d (%.1f%%)\n",
                length(mappable), length(mappable)/length(all_safe)*100))
    if (length(mappable) > 0) {
      cat("  Examples:\n")
      for (i in 1:min(5, length(mappable))) {
        orig <- name_map$Original_Name[name_map$Safe_Name == mappable[i]]
        cat(sprintf("    %s → %s\n", mappable[i], orig))
      }
    }
    unmapped <- setdiff(all_safe, name_map$Safe_Name)
    if (length(unmapped) > 0) {
      cat(sprintf("  Unmapped: %d metabolites\n", length(unmapped)))
      cat("  (These may already be original names or not present in mapping.)\n")
    }
  }
}

cat(paste0("\n", strrep("=", 50), "\n"))
cat("🎉 NAME MAPPING COMPLETE\n")
cat(paste0(strrep("=", 50), "\n"))

cat("Mapped files are saved as *_Mapped.csv in the output/ directory.\n")
cat("Use these for final reporting and figures with human-readable metabolite names.\n")
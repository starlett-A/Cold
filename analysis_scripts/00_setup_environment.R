# =================================================================
# 📦 SETUP ENVIRONMENT AND LOAD PACKAGES
# Purpose: Initialize R environment and load all required packages
# Note: This script should be run first!
# =================================================================

cat("🔧 Setting up analysis environment...\n")

# Install pacman if missing
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# Define required packages
required_packages <- c(
  # Data manipulation
  "dplyr", "tidyr", "tibble", "readr", "tidyverse",
  
  # Visualization
  "ggplot2", "ggrepel", "patchwork", "viridis", "RColorBrewer", "scales", "cowplot",
  
  # Statistics and Normalization
  "mixOmics", "imputeLCMD", "randomForest", "effsize", "pROC", "nortest",
  
  # Utilities
  "magrittr", "stringr", "matrixStats", "gridExtra"
)

# Load all packages (install if missing)
p_load(char = required_packages)

# Set random seed for reproducibility
set.seed(2026)

# Print confirmation
cat("✅ Environment setup complete\n")
cat("Loaded packages:\n")
print(sessionInfo()$otherPkgs)

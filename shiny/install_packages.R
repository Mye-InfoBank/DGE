#!/usr/bin/env Rscript

# Function to install packages if not already installed
install_if_missing <- function(packages, repos = "https://cloud.r-project.org/") {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg, repos = repos, dependencies = TRUE)
    }
  }
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(packages) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing Bioconductor package", pkg, "...\n")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }
}

# Essential CRAN packages for the Shiny app
cran_packages <- c(
  "shiny",           # Core Shiny framework
  "data.table",      # Fast data manipulation
  "ggplot2",         # Plotting
  "DT",              # DataTables for R
  "bslib",           # Bootstrap themes for Shiny
  "shinycssloaders", # Loading spinners
  "ggrepel",         # Text labels for plots
  "base64enc",       # Base64 encoding for images
  "future"           # Parallel processing
)

# Essential Bioconductor packages
bioc_packages <- c(
  "msigdbr"          # MSigDB gene sets
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install_if_missing(cran_packages)

# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")
install_bioc_if_missing(bioc_packages)

cat("All packages installed successfully!\n")

# Print session info for debugging
cat("\nR session info:\n")
print(sessionInfo()) 
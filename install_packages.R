packages <- c(
  "vegan",          # Ecological diversity analysis
  "randomForest",   # Machine learning
  "ggplot2",        # Visualization
  "reshape2",       # Data manipulation
  "pheatmap",       # Heatmaps
  "RColorBrewer",   # Color palettes
  "dplyr",          # Data manipulation
  "tidyr",          # Data tidying
  "jsonlite"        # JSON summaries
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

if (!require("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
}


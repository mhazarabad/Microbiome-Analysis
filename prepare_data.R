suppressPackageStartupMessages({
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    }
    BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
  }
  library(phyloseq)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", repos = "http://cran.us.r-project.org")
  }
  library(jsonlite)
})

dir.create("data", showWarnings = FALSE)

kostic_url <- "https://raw.githubusercontent.com/joey711/shiny-phyloseq/master/data/kostic.RData"
kostic_rdata <- file.path("data", "kostic.RData")

if (!file.exists(kostic_rdata)) {
  download.file(kostic_url, destfile = kostic_rdata, mode = "wb", quiet = TRUE)
}

load(kostic_rdata)

if (!inherits(kostic, "phyloseq")) {
  stop("Expected a phyloseq object named 'kostic' in the downloaded file.")
}

ps <- kostic
meta <- as(sample_data(ps), "data.frame")
otu <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
  otu <- t(otu)
}

meta$DIAGNOSIS <- as.character(meta$DIAGNOSIS)
meta$DIAGNOSIS[meta$DIAGNOSIS == "None"] <- NA

keep_samples <- rownames(meta)[!is.na(meta$DIAGNOSIS) & meta$DIAGNOSIS %in% c("Healthy", "Tumor")]

meta <- meta[keep_samples, , drop = FALSE]
otu <- otu[keep_samples, , drop = FALSE]

if (!all(rownames(meta) == rownames(otu))) {
  stop("Sample ordering mismatch between metadata and OTU table")
}

condition <- ifelse(meta$DIAGNOSIS == "Tumor", "Cancer", "Healthy")

parse_numeric <- function(x) {
  as.numeric(suppressWarnings(as.character(x)))
}

metadata <- data.frame(
  SampleID = rownames(meta),
  Condition = factor(condition, levels = c("Healthy", "Cancer")),
  Age = parse_numeric(meta$AGE),
  Sex = as.character(meta$SEX),
  Country = as.character(meta$COUNTRY),
  Treatment = as.character(meta$TREATMENT),
  stringsAsFactors = FALSE
)

metadata$Sex[metadata$Sex %in% c("None", "")] <- NA
metadata$Country <- sub("^GAZ:", "", metadata$Country)
metadata$Treatment[metadata$Treatment == "unaffected mucosa"] <- "Mucosa"
metadata$Treatment[metadata$Treatment == "tumor"] <- "Tumor"

otu_df <- cbind(SampleID = rownames(otu), as.data.frame(otu, check.names = FALSE))

taxonomy <- as.data.frame(tax_table(ps))
taxonomy$OTU <- rownames(taxonomy)
taxonomy <- taxonomy[, c("OTU", setdiff(names(taxonomy), "OTU"))]

write.csv(taxonomy, "data/taxonomy_table.csv", row.names = FALSE, quote = FALSE)

write.csv(metadata, "data/metadata.csv", row.names = FALSE, quote = FALSE)
write.csv(otu_df, "data/otu_table.csv", row.names = FALSE, quote = FALSE)
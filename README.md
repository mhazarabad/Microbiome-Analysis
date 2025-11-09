# Kostic 2012 Colorectal Cancer Microbiome

I swapped out the toy 16S OTU table for the real colorectal carcinoma cohort published by Kostic *et al.* (2012, *Genome Research*). The script `prepare_data.R` now downloads the public `kostic.RData` phyloseq object from Joey McMurdie's `shiny-phyloseq` repo, extracts 185 mucosal biopsies (95 healthy mucosa vs 90 tumor-adjacent tissue), and writes clean CSVs plus the original taxonomy table. Everything downstream runs on those actual counts and metadata.

## Data Provenance

- **Study:** Fusobacterium enrichment in colorectal carcinoma (Kostic *et al.*, 2012)
- **Accession:** NCBI SRA SRP018175 (amplicon region V3-V5)
- **Download script:** `prepare_data.R` (saves `data/otu_table.csv`, `data/metadata.csv`, `data/taxonomy_table.csv`)
- **Cohort summary:** 185 biopsies (median age 72), 2,505 OTUs after import, 119 OTUs retained after prevalence filtering (≥10% samples & ≥10 counts)

## Pipeline

1. Total-sum scaling to relative abundance after dropping rare taxa
2. Alpha diversity (Shannon, Simpson) with Wilcoxon tests
3. Bray–Curtis + PCA for beta diversity and PERMANOVA (999 permutations)
4. Taxon-wise Wilcoxon tests with BH FDR correction for differential abundance
5. Random Forest classifier (500 trees, 10 variables per split) with feature importance
6. Figure generation in `results/` + machine-readable outputs (`differential_abundance.csv`, `random_forest_importance.csv`, `rf_performance.txt`)

## Methods

**Normalization**: Total-sum scaling (TSS) to relative abundances

**Alpha Diversity**: Shannon and Simpson indices

**Beta Diversity**: Bray-Curtis dissimilarity with PCA ordination

**Differential Abundance**: Wilcoxon rank-sum test with Benjamini-Hochberg FDR correction

**Classification**: Random Forest with 10-fold cross-validation


## Results Snapshot

**Diversity shifts (Healthy vs Cancer)**

- Shannon p-value `4.44e-05`, Simpson p-value `9.09e-06` → tumors show significantly lower diversity; see `results/alpha_diversity.png` (box + jitter facets).

**Community structure**

- Bray–Curtis PERMANOVA p-value `0.001` (R² ≈ 0.05). PC1 separates tumor biopsies with Fusobacterium enrichment; see `results/beta_diversity_pca.png` for the PCA scatter with 95% ellipses.

**Differential abundance (FDR < 0.05)**

| OTU ID | Genus | log2FC (Cancer/Healthy) | Notes |
| --- | --- | --- | --- |
| 64396 | *Fusobacterium* | +1.92 | hallmark Kostic signal, strongly enriched in tumors |
| 374052 | *Fusobacterium* | +1.38 | second Fusobacterium OTU with elevated abundance |
| 72853 | *Faecalibacterium* | –1.74 | anti-inflammatory commensal depleted in tumors |
| 322235 | *Bacteroides* | –1.47 | healthy mucosa enriched |

The full ranked table is in `results/differential_abundance.csv`, and I plotted the top 15 FDR-significant features in `results/differential_abundance_heatmap.png`.

**Random Forest sanity check**

- OOB accuracy `55.1%` (error 44.9%) → the simple abundance signature alone is noisy (matching the Kostic paper’s note that taxonomy-only models are imperfect).
- Top importance scores (`results/random_forest_importance.csv` / `results/random_forest_importance.png`) still elevate the same *Fusobacterium* OTUs, even though the model struggles to predict subject-level labels.

## Literature cross-check

- **Fusobacterium enrichment:** Kostic et al. reported striking Fusobacterium over-representation in tumors. My Wilcoxon analysis recovers multiple Fusobacterium OTUs with log₂ fold-changes > +1.3 (FDR < 0.04) and the PCA plot shows the same axis of separation.
- **Loss of butyrate producers:** Healthy mucosa is enriched for *Faecalibacterium*, *Ruminococcus*, and *Collinsella*—all depleted in tumors with FDR-adjusted p-values < 0.001—reproducing the butyrate depletion motif from the paper.
- **Predictive modeling:** The Kostic study combined microbial features with host data for better discrimination; my pure-OTU Random Forest mirrors their observation that taxonomy alone is insufficient (OOB accuracy barely >50%).

## How to run the project

```bash
Rscript install_packages.R
Rscript prepare_data.R
Rscript microbiome_analysis.R
```

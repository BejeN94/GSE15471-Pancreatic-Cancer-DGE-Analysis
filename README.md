# GSE15471 pancreatic cancer differential expression

Reproducible analysis of the public GEO microarray series [GSE15471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471), which profiles pancreatic ductal adenocarcinoma (PDAC) and non-tumor pancreatic tissue on the Affymetrix Human Genome U133 Plus 2.0 array (GPL570).

The GEO series reports paired tumor and normal specimens. The workflow downloads the processed GEO matrix, writes a sample annotation template for review, and then uses limma to compare tumor with normal tissue. It does not treat this microarray as RNA-seq count data.

## Requirements

- R (4.2 or newer recommended)
- Internet access on the first run
- Bioconductor packages `GEOquery` and `limma`

Install the packages once:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"))
```

## Run

From the repository root:

```r
source("analysis.R")
```

On the first run, the script downloads GSE15471 into `data/raw/`, exports the GEO sample metadata to `results/sample_annotations.csv`, and stops. Review that CSV before analysis:

- Set `group` to exactly `tumor` or `normal` for each included sample.
- Set `include` to `TRUE` or `FALSE` to exclude a sample.
- Fill `pair_id` for matched tumor/normal samples using the same patient identifier for both samples. Leave it blank for every sample if reliable patient IDs are unavailable.
- Do not infer pairing from row order. Confirm labels and pair IDs against GEO sample records.

Run `source("analysis.R")` again after saving the reviewed CSV. The script validates the annotations and fails with a specific message when groups or pairs are incomplete.

## Outputs

Files are written to `results/`:

- `sample_annotations.csv`: editable sample labels and optional matched-pair IDs.
- `sample_metadata.csv`: GEO sample information used in the run.
- `expression_matrix.rds`: downloaded processed probe-level matrix.
- `differential_expression_probes.csv`: all probe-level limma results, including log2 fold change, average expression, moderated t statistic, raw p-value and Benjamini–Hochberg adjusted p-value.
- `differential_expression_genes.csv`: probe results with gene symbols where GEO platform annotation provides them.
- `volcano_plot.png`: overview of probe-level differential expression.
- `sessionInfo.txt`: R and package versions for reproducibility.

Positive log2 fold change means higher expression in tumor. Significance is reported using adjusted p-values; the default exploratory threshold is FDR < 0.05. Probe-level findings should be interpreted with the array platform and tissue composition in mind.

## Methods and provenance

Expression values come from the processed series matrix supplied by GEO, rather than raw CEL files. If the downloaded matrix is not already on a log2 scale, the script applies log2 transformation after checking that values are positive. The linear model uses limma empirical Bayes moderation. If complete pair IDs are supplied, the model includes patient as a blocking factor; otherwise it fits an unpaired tumor-versus-normal comparison.

The GEO record describes the cohort and platform: [GSE15471 at NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471). This repository provides an analysis workflow, not clinical or diagnostic guidance.

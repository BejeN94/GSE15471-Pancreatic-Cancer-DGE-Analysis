# Differential expression analysis for GEO series GSE15471
# Run from the repository root with: source("analysis.R")

options(stringsAsFactors = FALSE)

required_packages <- c("GEOquery", "limma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install missing Bioconductor package(s) first: ",
    paste(missing_packages, collapse = ", "),
    ". See the README for installation instructions.",
    call. = FALSE
  )
}

dir.create(file.path("data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

geo_id <- "GSE15471"
annotation_file <- file.path("results", "sample_annotations.csv")

# GEOquery caches downloaded files in GSEMatrix/soft format under data/raw.
geo <- GEOquery::getGEO(
  geo_id,
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = file.path("data", "raw")
)
if (!length(geo)) stop("GEOquery returned no expression sets for ", geo_id, ".")
eset <- geo[[1]]
expression <- Biobase::exprs(eset)
pheno <- Biobase::pData(eset)

if (is.null(colnames(expression)) || is.null(rownames(expression))) {
  stop("The GEO expression matrix is missing sample or probe identifiers.")
}
if (!all(colnames(expression) %in% rownames(pheno))) {
  stop("Expression matrix samples do not all match GEO phenotype records.")
}
pheno <- pheno[colnames(expression), , drop = FALSE]

collapse_text <- function(data, columns) {
  columns <- intersect(columns, names(data))
  if (!length(columns)) return(rep("", nrow(data)))
  apply(data[, columns, drop = FALSE], 1, function(row) {
    paste(as.character(row[!is.na(row) & nzchar(as.character(row))]), collapse = " | ")
  })
}

char_cols <- grep("^characteristics_ch1", names(pheno), value = TRUE)
label_text <- tolower(paste(
  pheno$title,
  if ("source_name_ch1" %in% names(pheno)) pheno$source_name_ch1 else "",
  collapse_text(pheno, char_cols)
))
guess_group <- ifelse(
  grepl("normal|non[- ]?tumou?r|control|healthy|adjacent", label_text), "normal",
  ifelse(grepl("tumou?r|pdac|carcinoma|cancer", label_text), "tumor", "")
)

source_name <- if ("source_name_ch1" %in% names(pheno)) {
  as.character(pheno$source_name_ch1)
} else {
  rep("", nrow(pheno))
}

# A first run creates a reviewable template, never silently proceeds on guesses.
if (!file.exists(annotation_file)) {
  annotations <- data.frame(
    geo_accession = rownames(pheno),
    title = as.character(pheno$title),
    source_name = source_name,
    group = guess_group,
    pair_id = "",
    include = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.csv(annotations, annotation_file, row.names = FALSE, na = "")
  utils::write.csv(pheno, file.path("results", "sample_metadata.csv"), row.names = FALSE)
  saveRDS(expression, file.path("results", "expression_matrix.rds"))
  stop(
    "Created results/sample_annotations.csv. Review group, include, and pair_id ",
    "against GEO records, save the file, then run source('analysis.R') again.",
    call. = FALSE
  )
}

annotations <- utils::read.csv(annotation_file, colClasses = "character",
                               na.strings = character(), check.names = FALSE)
required_columns <- c("geo_accession", "group", "pair_id", "include")
missing_columns <- setdiff(required_columns, names(annotations))
if (length(missing_columns)) {
  stop("Sample annotation file is missing: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(annotations$geo_accession)) {
  stop("geo_accession values in sample_annotations.csv must be unique.")
}
if (!setequal(annotations$geo_accession, colnames(expression))) {
  missing <- setdiff(colnames(expression), annotations$geo_accession)
  extra <- setdiff(annotations$geo_accession, colnames(expression))
  stop(
    "Annotations do not match expression samples. Missing: ",
    paste(missing, collapse = ", "), "; extra: ", paste(extra, collapse = ", "),
    ". Delete the template and rerun to regenerate it."
  )
}
annotations <- annotations[match(colnames(expression), annotations$geo_accession), ]
included <- tolower(trimws(annotations$include)) %in% c("true", "t", "1", "yes", "y")
if (!any(included)) stop("No samples are marked include=TRUE.")
annotations <- annotations[included, , drop = FALSE]
expression <- expression[, included, drop = FALSE]

group <- tolower(trimws(annotations$group))
if (any(!group %in% c("tumor", "normal"))) {
  stop("Every included sample must have group exactly 'tumor' or 'normal'.")
}
group <- factor(group, levels = c("normal", "tumor"))
if (any(table(group) < 2L)) {
  stop("At least two included samples per group are required.")
}

pair_id <- trimws(annotations$pair_id)
has_pairs <- any(nzchar(pair_id))
if (has_pairs) {
  if (any(!nzchar(pair_id))) {
    stop("For paired analysis, fill pair_id for every included sample or leave all blank.")
  }
  pair_id <- factor(pair_id)
  pair_group_counts <- table(pair_id, group)
  if (any(pair_group_counts != 1L)) {
    stop("Each pair_id must have exactly one normal and one tumor sample.")
  }
  design <- stats::model.matrix(~ pair_id + group)
  contrast_name <- "grouptumor"
} else {
  design <- stats::model.matrix(~ group)
  contrast_name <- "grouptumor"
}
if (qr(design)$rank < ncol(design)) stop("The design matrix is not full rank.")

# GEO processed matrices are often already log2 transformed. Transform only
# clearly positive, unlogged intensity-scale values.
if (any(!is.finite(expression))) stop("Expression matrix contains non-finite values.")
if (min(expression, na.rm = TRUE) > 0 &&
    max(expression, na.rm = TRUE) > 100) {
  expression <- log2(expression)
}
if (min(expression, na.rm = TRUE) < 0 && max(expression, na.rm = TRUE) > 100) {
  stop("Expression scale is ambiguous; inspect the GEO matrix before analysis.")
}

fit <- limma::lmFit(expression, design)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
if (!contrast_name %in% colnames(fit$coefficients)) {
  stop("Could not find tumor-versus-normal coefficient in the fitted model.")
}
probe_results <- limma::topTable(
  fit,
  coef = contrast_name,
  number = Inf,
  sort.by = "P",
  adjust.method = "BH"
)
probe_results$probe_id <- rownames(probe_results)
probe_results <- probe_results[, c(
  "probe_id", setdiff(names(probe_results), "probe_id")
), drop = FALSE]
utils::write.csv(
  probe_results,
  file.path("results", "differential_expression_probes.csv"),
  row.names = FALSE
)

# Attach available symbols from the GEO platform feature annotation.
feature_data <- Biobase::fData(eset)
symbol_col <- grep("gene.?symbol|symbol", names(feature_data), ignore.case = TRUE, value = TRUE)
if (length(symbol_col) && !is.null(rownames(feature_data))) {
  feature_idx <- match(probe_results$probe_id, rownames(feature_data))
  symbols <- rep(NA_character_, nrow(probe_results))
  has_feature <- !is.na(feature_idx)
  symbols[has_feature] <- as.character(
    feature_data[[symbol_col[1]]][feature_idx[has_feature]]
  )
  probe_results$gene_symbol <- symbols
} else {
  probe_results$gene_symbol <- NA_character_
}
utils::write.csv(
  probe_results,
  file.path("results", "differential_expression_genes.csv"),
  row.names = FALSE
)

plot_data <- probe_results
plot_data$neg_log10_fdr <- -log10(pmax(plot_data$adj.P.Val, .Machine$double.xmin))
plot_data$significant <- plot_data$adj.P.Val < 0.05
grDevices::png(file.path("results", "volcano_plot.png"), width = 1600, height = 1200, res = 180)
graphics::plot(
  plot_data$logFC, plot_data$neg_log10_fdr,
  pch = 16, cex = 0.45,
  col = ifelse(plot_data$significant, "#C43C39", "#666666"),
  xlab = "log2 fold change (tumor vs normal)",
  ylab = "-log10 Benjamini-Hochberg adjusted p-value",
  main = paste(geo_id, "differential expression")
)
graphics::abline(h = -log10(0.05), lty = 2, col = "#444444")
grDevices::dev.off()

utils::write.csv(pheno, file.path("results", "sample_metadata.csv"), row.names = FALSE)
saveRDS(expression, file.path("results", "expression_matrix.rds"))
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))

message(
  "Analysis complete: ", sum(group == "tumor"), " tumor and ",
  sum(group == "normal"), " normal samples; ",
  if (has_pairs) length(levels(pair_id)) else "unpaired",
  " model. Results are in results/."
)

# log <- file(snakemake@log[[1]], open="wt")
# sink(log)
# sink(log, type="message")

library("DESeq2")
parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}
if (!exists("snakemake")) {
  library(magrittr)
  BASE_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/apelin_exp"
  cts <- file.path(BASE_DIR, "counts/all.tsv") %>%
    read.table(header = T, row.names = "gene", check.names = F)
  coldata <- read.table("../data/apelin_2020/VascAge_samples_apelin.tsv",
    header = TRUE,
    row.names = "sample",
    check.names = FALSE
  )
}
# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]],
  header = TRUE,
  row.names = "gene",
  check.names = FALSE
)
coldata <- read.table(snakemake@params[["samples"]],
  header = TRUE,
  row.names = "sample", check.names = FALSE
)



## Reorder coldata rows to match cts col order (beacause deseq things)
if (!all(colnames(cts) == rownames(coldata))) {
  sample_ids <- colnames(cts)
  coldata <- coldata[sample_ids, ]
}

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  # design = as.formula("~condition"))
  design = as.formula(snakemake@params[["model"]]),
)




# remove genes defined in config
excluded_genes <- snakemake@config[["excluded_genes"]]
# print(excluded_genes)
if (length(excluded_genes) > 0) {
  ex_index <- purrr::map_int(excluded_genes,
    grep, rownames(dds),
    value = F
  )
  if (length(ex_index) > 0) {
    dds <- dds[-ex_index, ]
  }
}
# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 9, ]
# normalization and preprocessing
dds <- DESeq(dds,
  parallel = parallel
)

saveRDS(dds, file = snakemake@output[[1]])

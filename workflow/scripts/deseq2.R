log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

contrast <- c(snakemake@wildcards[["condition"]], snakemake@params[["contrast"]])
if (!is.null(snakemake@config$diffexp$custom_model[[contrast[1]]])) {
  # Rerun deseq since it was intialized 
  colData(dds)[,contrast[1]] <- as.factor(colData(dds)[,contrast[1]])
  design(dds) <- as.formula(snakemake@config$diffexp$custom_model[[contrast[1]]])
  dds <- DESeq(dds)
}
res <- results(dds, contrast = contrast, parallel = parallel)
resstat <- res$stat
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")
# sort by p-value
res$stat <- resstat
res <- res[, c(
  "baseMean",
  "log2FoldChange",
  "lfcSE", "stat",
  "pvalue", "padj"
)]
res <- res[order(res$padj), ]

# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim = c(-2, 2))
dev.off()

write.table(as.data.frame(res),
  sep = "\t",
  quote = FALSE, file = snakemake@output[["table"]]
)

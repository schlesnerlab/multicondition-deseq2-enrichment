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

lfc_shrink <- snakemake@config[["diffexp"]][["shrink_lfc"]]
dds <- readRDS(snakemake@input[[1]])
cond_id <- snakemake@wildcards[["condition"]]
contrast_groups <- snakemake@params[["contrast"]]
contrast <- c(cond_id, contrast_groups)
print(contrast)

coef_string <- paste0(cond_id, "_", contrast_groups[[1]], "_vs_", contrast_groups[[2]])
print(coef_string)
if (!is.null(snakemake@config$diffexp$custom_model[[contrast[1]]])) {
  # Rerun deseq since it was intialized
  for (x in names(snakemake@config[["diffexp"]][["contrasts"]])) {
    colData(dds)[, x] <- as.factor(colData(dds)[, x])
  }

  colData(dds)[, contrast[1]] <- as.factor(colData(dds)[, contrast[1]])
  design(dds) <- as.formula(snakemake@config$diffexp$custom_model[[contrast[1]]])
  dds <- DESeq(dds)
}

zero_index <- apply(counts(dds)[, colData(dds)[, cond_id] %in% contrast_groups], 1, function(x) {
  length(which(x > 0))
})
# dds <- dds[zero_index > ncol(dds[,colData(dds)[,cond_id] %in% contrast_groups])/2,]
res <- results(dds, contrast = contrast, parallel = parallel)
resstat <- res$stat
# shrink fold changes for lowly expressed genes
if (lfc_shrink) {
  dds@colData[, cond_id] <- relevel(dds@colData[, cond_id], ref = contrast_groups[2])
  dds <- DESeq(dds, parallel = parallel)
  res_new <- lfcShrink(dds, coef = stringr::str_replace_all(coef_string, "-", "."), type = "apeglm")
  print(head(res_new))
  # sort by p-value
  res_new$stat <- resstat
  res_new <- res_new[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  res_new <- res_new[order(res_new$padj), ]
}
# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim = c(-2, 2))
if (lfc_shrink) {
  plotMA(res_new, ylim = c(-2, 2))
}
dev.off()

if (lfc_shrink) {
  res <- res_new
}
#


write.table(as.data.frame(res), sep = "\t", quote = FALSE, file = snakemake@output[["table"]])

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
print(contrast)
if (snakemake@config[["perform_perturbation"]]) {
  # Res
  index <- (dds@colData$condition %in% snakemake@params[["contrast"]])
  cond_labels <- dds@colData$condition[index]
  pert_cond_labels <- cond_labels
  while (isTRUE(all.equal(pert_cond_labels, cond_labels))) {
    pert_cond_labels <- sample(cond_labels,
      length(cond_labels),
      replace = F
    )
  }
  print(pert_cond_labels)
  dds@colData[index, "condition"] <- pert_cond_labels
}
#dds$condition <- relevel(dds$condition, contrast[3])
res <- results(dds, contrast = contrast, parallel = parallel)
resstat <- res$stat
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast = contrast, res = res,type = "ashr")
# sort by p-value
res$stat <- resstat
res <- res[,c(
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

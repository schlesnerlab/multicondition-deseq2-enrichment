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

res <- results(dds, contrast = contrast, parallel = parallel)
# shrink fold changes for lowly expressed genes
# LFC shrinkage remains deactivated since it doesn't work with support how we run
# DESeq2
# 
#res <- lfcShrink(dds, contrast = contrast , res = res)

# sort by p-value
res <- res[order(res$padj), ]

# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim = c(-2, 2))
dev.off()

write.table(as.data.frame(res),
  sep = "\t",
  quote = FALSE, file = snakemake@output[["table"]]
)

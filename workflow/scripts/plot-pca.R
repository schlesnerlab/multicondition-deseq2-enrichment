log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
#counts <- assay(dds)
svg(snakemake@output[[1]])
plotPCA(dds, intgroup = snakemake@params[["pca_labels"]])
dev.off()

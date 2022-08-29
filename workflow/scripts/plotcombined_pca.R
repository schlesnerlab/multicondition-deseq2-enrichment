library("DESeq2")
library(ReportingTools)
library(tidyverse)
# load deseq2 data
BASE_ANALYSIS_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/"
vage_dds_path <- file.path(paste0(BASE_ANALYSIS_DIR), "Vasc_age2020/deseq2/all.rds")
ape_dds_path <- file.path(BASE_ANALYSIS_DIR, "apelin_exp/deseq2/all.rds")

vag_dds <- readRDS(vage_dds_path)
ape_dds <- readRDS(ape_dds_path)
fpkm_path <- file.path(BASE_ANALYSIS_DIR, "Vasc_age2020/fpkm/all.tsv")

fpkm_t <- readr::read_tsv(fpkm_path)

# Get young + aged groups

vag_coldata <- colData(vag_dds) %>%
  as_tibble(rownames = "id") %>%
  add_column(batch = "vascage2020")
ape_coldata <- colData(ape_dds) %>%
  as_tibble(rownames = "id") %>%
  add_column(batch = "apelin_exp")

all_coldata <- rbind(vag_coldata, ape_coldata)
my_col <- all_coldata %>%
  dplyr::filter(condition %in% c("d0-lung", "18m-lung", "Young-EC", "Aged-EC")) %>%
  dplyr::arrange(batch, condition)

my_col[my_col$condition == "18m-lung", "condition"] <- "Aged-EC"
my_col[my_col$condition == "d0-lung", "condition"] <- "Young-EC"
## Counts

vag_c <- counts(vag_dds) %>% as_tibble(rownames = "geneid")
ape_c <- counts(ape_dds) %>% as_tibble(rownames = "geneid")

merged_counts <- dplyr::inner_join(vag_c, ape_c, )
gene_ids <- merged_counts$geneid
count_mat <- as.matrix(merged_counts[, -c(1)], )

rownames(count_mat) <- gene_ids
count_mat <- count_mat[, my_col$id]

dds <- DESeqDataSetFromMatrix(count_mat, colData = my_col, design = ~ condition + batch)



# obtain normalized counts
counts <- rlog(dds, blind = FALSE)

plotPCA(counts, intgroup = c("batch", "condition")) + theme_bw()

library(PCAtools)
yet <- PCAtools::pca(assay(counts), metadata = colData(counts))

screeplot(yet)

biplot(yet, colby = "condition", shape = "batch", legendPosition = "right", lab = NULL, drawConnectors = F)

dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
res <- DESeq2::results(dds, contrast = c("condition", "Young-EC", "Aged-EC"), tidy = T) %>% dplyr::rename(c("gene_id" = "row"))
filer <- fpkm_t %>% dplyr::filter(gene %in% res$gene_id)
library(RNAscripts)
joined_df <- join_tables(res, filer)
res_orig <- DESeq2::results(dds, contrast = c("condition", "Young-EC", "Aged-EC"))
EnhancedVolcano::EnhancedVolcano(as.data.frame(joined_df),
  lab = joined_df$gname,
  x = "log2FoldChange",
  y = "padj",
  ylab = bquote(~ -Log[10] ~ italic(Padj)), FCcutoff = 1
)

plotMA(res_orig)

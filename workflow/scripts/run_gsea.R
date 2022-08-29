library(clusterProfiler)


## Snakemake header

if (exists("snakemake")) {
  diffexp_tb_path <- snakemake@input[["table"]]
  fpkm_path <- snakemake@input[["fpkm_path"]]
  contrast_groups <- snakemake@params[["contrast"]]
  pvalue_threshold <- snakemake@config[["diffexp"]][["pval_threshold"]]
  LFC_threshold <- snakemake@config[["diffexp"]][["LFC_threshold"]]
  out_file <- snakemake@output[["gsea_result"]]
} else {
  BASE_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/apelin_exp"
  contrast <- "AEC-vs-YEC"
  contrast_groups <- c("AEC", "YEC")
  diffexp_tb_path <- file.path(
    BASE_DIR,
    glue::glue("results/diffexp/{contrast}.diffexp.tsv")
  )
  fpkm_path <- file.path(BASE_DIR, "fpkm/all.tsv")
  pvalue_threshold <- 0.05
  LFC_threshold <- 0.5
}


diffxp_tb <- readr::read_tsv(diffexp_tb_path,
  col_names = c(
    "gene_id",
    "baseMean",
    "logFoldChange",
    "lfcSE", "stat",
    "pvalue", "padj"
  ),
  skip = 1
)
fpkm <- readr::read_tsv(fpkm_path)

# Read Files
filer <- fpkm %>% dplyr::filter(gene %in% diffxp_tb$gene_id)

joined_df <- RNAscripts::join_tables(diffxp_tb, filer)

joined_df <- joined_df %>%
  dplyr::mutate(overexpressed_in = ifelse(logFoldChange > 0,
    contrast_groups[1],
    contrast_groups[2]
  ))

joined_df$gsea_stat <- joined_df$stat
gene_list <- joined_df %>% dplyr::select(c(gname, gsea_stat))
ensemblgene_list <- joined_df %>% dplyr::select(c(gene, gsea_stat))
de_genes <- joined_df %>%
  dplyr::filter(padj < pvalue_threshold & abs(gsea_stat) > LFC_threshold) %>%
  dplyr::select(c(gene, stat, gsea_stat))


t_table <- RNAscripts::get_entrezgenes_from_ensembl(filer %>% dplyr::pull(gene),
  input_type = "ENSEMBL"
) %>%
  RNAscripts::table_to_list(., "ENTREZID", "ENSEMBL")

msig_enrichment <- RNAscripts::run_msig_enricher(list(de_genes),
  universe = ensemblgene_list$gene, GSEA = FALSE,
  translation_table = t_table, msdb_var = "entrez_gene",
  input_type = "ENSEMBL", category = "H"
)[[1]]
# Run the MSIG enricher
msig_gsea <- RNAscripts::run_msig_enricher(list(ensemblgene_list),
  translation_table = t_table,
  msdb_var = "entrez_gene",
  input_type = "ENSEMBL",
  category = "H", eps = 0
)[[1]]

msig_c3 <- RNAscripts::run_msig_enricher(list(ensemblgene_list),
  translation_table = t_table,
  subcategory = "TFT:GTRD",
  msdb_var = "entrez_gene",
  input_type = "ENSEMBL",
  category = "C3", eps = 0
)[[1]]

msig_c6 <- RNAscripts::run_msig_enricher(list(ensemblgene_list),
  translation_table = t_table,
  msdb_var = "entrez_gene",
  input_type = "ENSEMBL",
  category = "C6", eps = 0
)[[1]]
kegg <- RNAscripts::run_gsea(ensemblgene_list,
  input_type = "ENSEMBL", p_valcut = 0.1
)

g_vec <- RNAscripts::get_entrezgene_vector(ensemblgene_list, "ENSEMBL")
reactome_stuff <- ReactomePA::gsePathway(g_vec,
  organism = "mouse",
  verbose = TRUE
)

enrich_list <- list(
  "msig_enrichment" = msig_enrichment,
  "msig_gsea" = msig_gsea,
  "msig_C3" = msig_c3,
  "msig_C6" = msig_c6,
  "kegg" = kegg,
  "reactome_stuff" = reactome_stuff
)
saveRDS(
  object = enrich_list,
  file = out_file
)

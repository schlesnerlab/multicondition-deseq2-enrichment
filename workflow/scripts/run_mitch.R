# This script runs mitch MANOVA gene set enrichment across multiple DEseq2 runs Writes the results to a tsv file and
# renders a report for the run experiment.  TOdo setup which pipelines this should run with

library(tibble)
if (!require(RNAscripts)) {
  devtools::install("./scripts/RNAscripts", upgrade = "never")
}
library(RNAscripts)
if (packageVersion("mitch") < "1.1.8") {
  devtools::install_github("markziemann/mitch")
}
# library(org.Mm.eg.db)
library(mitch)
library(tidyverse)
library(ensembldb)
library(AnnotationHub)
## Load the annotation resource.
ah <- AnnotationHub()
Gtf <- ah["AH28753"]
## Create a EnsDb database file from this.
DbFile <- ensDbFromAH(Gtf)
## We can either generate a database package, or directly load the data
edb <- EnsDb(DbFile)
library(msigdbr)
if (exists("snakemake")) {
  # Input Files
  diffexp_tables_paths <- snakemake@input[["tables"]]
  contrast_list <- snakemake@params[["contrasts"]]
  fpkm_path <- snakemake@input[["fpkm"]]
  # Output Files
  mitch_table <- snakemake@output[["mitch_table"]]
  mitch_rds <- snakemake@output[["mitch_rds"]]
  report_path <- snakemake@output[["mitch_report"]]
  # Parameters
  threads <- snakemake@threads
  enrichment_term_type <- "SYMBOL"
} else {
  configfile <- yaml::read_yaml("./configs/VascAge_Apelin_config_wo_youngplus.yaml")
  BASE_ANALYSIS_DIR <- configfile$dirs$BASE_ANALYSIS_DIR
  contrast_list <- names(configfile$diffexp$contrasts$condition)
#  contrast_list <- c(
#    "basal_cre_pos_vs_basal_cre_neg", "tumor_cre_pos_vs_tumor_cre_neg", "tumor_cre_pos_vs_basal_cre_pos",
#    "tumor_cre_neg_vs_basal_cre_neg"
#  )
  diffexp_tables_paths <- as.list(file.path(BASE_ANALYSIS_DIR, glue::glue("results/diffexp/condition/{contrast_list}.diffexp.tsv")))


  fpkm_path <- file.path(BASE_ANALYSIS_DIR, "fpkm/all.tsv")
  threads <- 1
  report_path <- "big_yeeter.html"
  enrichment_term_type <- "SYMBOL"
}

## Read data
diff_exp_tables <- purrr::map(diffexp_tables_paths, readr::read_tsv, col_names = c(
  "gene_id", "baseMean", "logFoldChange",
  "lfcSE", "stat", "pvalue", "padj"
), skip = 1)
print(contrast_list)
names(diff_exp_tables) <- names(contrast_list)
diff_exp_frames <- purrr::map(diff_exp_tables, function(x) {
  gene_ids <- x$gene_id
  diff_exp_fr <- as.data.frame(x[, -1])
  rownames(diff_exp_fr) <- gene_ids
  diff_exp_fr
})
# names(diff_exp_frames) <- contrast_list Read FPKM file
fpkm <- readr::read_tsv(fpkm_path)

match_table <- fpkm %>%
  dplyr::select(c("gene", "gname"))
match_vec <- setNames(object = match_table$gname, nm = match_table$gene)


### Setupt mitch input data
prep_mitch_input <- function(deseq_list, e_term = "ENSEMBL") {
  # mitch_input_df <- mitch::mitch_import(diff_exp_frames, 'DESeq2', geneTable = as.data.frame(fpkm[,c('gene',
  # 'gname')]))
  mitch_input_df <- mitch::mitch_import(deseq_list, "DESeq2")
  gname_tb <- tibble::tibble(Transcript_id = rownames(mitch_input_df))
  if (e_term == "ENSEMBL") {
    gname_tb["gene_id"] <- stringr::str_extract(gname_tb$Transcript_id, "ENSMUSG[0-9]*")
  } else if (e_term == "SYMBOL") {
    gname_tb["gene_id"] <- match_vec[gname_tb$Transcript_id]
  }
  rownames(mitch_input_df) <- gname_tb$gene_id
  colnames(mitch_input_df) <- stringr::str_remove_all(colnames(mitch_input_df), pattern = "-")
  ifelse(anyDuplicated(gname_tb$gene_id), warning("Duplicated genes in mitch input dataframe."), print("no duplicates"))

  mitch_input_df
}
#' Uniquify a vector
#'
#' @param x Input factor to check
#' @return corrected facotr value
#' @examples
#' NULL
uniquify <- function(x) {
  while (anyDuplicated(x)) {
    x[duplicated(x)] <- paste0(x[duplicated(x)], "_1")
  }
  x
}
mitch_input_df <- prep_mitch_input(diff_exp_frames)
# mm_reactome <- buildReactome(output_type = enrichment_term_type) ensembl_reactome <- buildReactome(output_type =
# 'ENSEMBL')


if (enrichment_term_type == "SYMBOL") {
  # Get GENE ids
  new_ids <- ensembldb::select(edb, keys = rownames(mitch_input_df), keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
  update_vector <- setNames(new_ids$SYMBOL, nm = new_ids$GENEID)
  # Stop if any duplicated genes in reactome gene set
  names(match_vec) <- stringr::str_extract(names(match_vec), "ENSMUSG[0-9]*")

  geneIDs <- match_vec[rownames(mitch_input_df)]
  geneIDs[names(update_vector)] <- update_vector


  col_names <- c("gs_name", "gene_symbol")
  mm_msigdbr <- msigdbr::msigdbr(species = "Mus musculus", category = "H") %>%
    dplyr::select(col_names)
  msigdb_gene_list <- split(x = mm_msigdbr$gene_symbol, f = mm_msigdbr$gs_name)
  ensembl_in_msigdb <- new_ids %>%
    dplyr::filter(SYMBOL %in% mm_msigdbr$gene_symbol) %>%
    pull(GENEID)

  ## Cleanup duplicates
  mitch_input_df <- mitch_input_df[(rownames(mitch_input_df) %in% ensembl_in_msigdb) | !(duplicated(geneIDs[rownames(mitch_input_df)]) |
    duplicated(geneIDs[rownames(mitch_input_df)], fromLast = T)), ]

  geneIDs <- match_vec[rownames(mitch_input_df)]
  # geneIDs[names(update_vector)] <- update_vector
  dup_mitch <- mitch_input_df[] %>%
    as_tibble(rownames = "ENSEMBL")
  dup_mitch["SYMBOL"] <- geneIDs[dup_mitch$ENSEMBL]
  # Remove samples which(dup_mitch$ENSEMBL %in% unlist(ensembl_reactome))
  final_name_vec <- setNames(dup_mitch$SYMBOL, nm = dup_mitch$ENSEMBL)
  final_name_vec <- final_name_vec[-which(!(final_name_vec %in% mm_msigdbr$gene_symbol) & duplicated(final_name_vec))]
  # final_name_vec['ENSMUSG00000051396'] <- 'Gm45902'
  if (anyDuplicated(final_name_vec)) {
    stop("Duplicate in gene names")
  }
  # geneIDs <- uniquify(geneIDs)
  mitch_input_df <- mitch_input_df[names(final_name_vec), ]
  rownames(mitch_input_df) <- final_name_vec
}

## Setup Databases to test msig_hallmarck <- msigdbr::msigdbr(species = 'Mus musculus')

# bain <- msig_hallmarck %>% dplyr::select(c('gs_name', 'gene_symbol')) %>% group_split(gs_name)

# more_bain <- purrr::map(bain, function(x){x$gene_symbol})

# mitch::mitch_calc(mitch_input_df, more_bain, cores = 6, priority ='significance')

# mm_reactome <- buildReactome(output_type = enrichment_term_type) print(colnames(mitch_input_df))
# print(anyDuplicated(colnames(mitch_input_df)))
reac_test <- mitch::mitch_calc(mitch_input_df, msigdb_gene_list, cores = threads)

readr::write_tsv(reac_test$enrichment_result, mitch_table)
saveRDS(reac_test, file = mitch_rds)

mitch::mitch_report(reac_test, outfile = report_path)

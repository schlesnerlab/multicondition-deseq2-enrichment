library(openxlsx)
library(magrittr)
library(dplyr)
library(tidyr)
library(glue)
if (!require(RNAscripts)) {
  devtools::install("./scripts/RNAscripts", upgrade = "never")
}
library(RNAscripts)
if (exists("snakemake")) {
  diffexp_tb_path <- snakemake@input[["table"]]
  cond_id <- snakemake@wildcards[["condition"]]
  fpkm_path <- snakemake@input[["fpkm"]]
  samp_map <- snakemake@params[["samp_map"]]
  contrast_groups <- snakemake@config[["diffexp"]][["contrasts"]][[cond_id]]
  contrast_names <- snakemake@params[["contrast_groups"]]
  names(contrast_groups) <- names(contrast_names)
  output_path <- snakemake@output[["outpath"]]
  print(c(
    diffexp_tb_path, fpkm_path, samp_map,
    contrast_groups, contrast_names
  ))
  pval_threshold <- snakemake@config[["diffexp"]][["pval_threshold"]]
  lfc_threshold <- snakemake@config[["diffexp"]][["LFC_threshold"]]
} else {
  BASE_ANALYSIS_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/apelin_exp_test"

  test_config <- yaml::read_yaml("/desktop-home/heyer/projects/Vascular_Aging/RNAseq/multicondition-deseq2-enrichment/configs/VascAge_Apelin_config.yaml")
  cond_id <- "condition"
  contrast_groups <- test_config$diffexp$contrasts$condition
  contrast_names <- test_config$diffexp$contrasts$condition
  diffexp_tb_path <- as.list(file.path(
    BASE_ANALYSIS_DIR,
    glue::glue("results/diffexp/condition/{names(contrast_groups)}.diffexp.tsv"
  )))

  names(diffexp_tb_path) <- names(contrast_groups)
  lfc_threshold <- 1
  pval_threshold <- 0.05

  fpkm_path <- file.path(BASE_ANALYSIS_DIR, "fpkm/all.tsv")
  samp_map <- "/desktop-home/heyer/projects/Vascular_Aging/RNAseq/rna-seq-star-deseq2/data/apelin_2020/VascAge_samples_apelin.tsv"
  output_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/apelin_exp_test/results/diffexp/condition/excel_tables/diff_exp_man.xlsx"
}


print("read samples")
sample_overview <- readr::read_tsv(samp_map)


col_name <- c("gene", "gname")

print("read fpkm table")
fpkm <- readr::read_tsv(fpkm_path)

build_output_table <- function(diff_exp_table, c_group, fpkm, col_name) {
  print(c_group)
  print(col_name)
  relevant_samples <- sample_overview %>%
    dplyr::filter(condition %in% c_group) %>%
    pull(sample)
  print(relevant_samples)
  deseq2_table <- readr::read_tsv(diff_exp_table,
    col_names = c(
      "gene_id",
      "baseMean",
      "logFoldChange",
      "lfcSE", "stat",
      "pvalue", "padj"
    ),
    skip = 1
  )
  print("read deseq")

  # print(colnames(fpkm))
  fpkm_data <- fpkm %>%
    dplyr::filter(gene %in% deseq2_table$gene_id) %>%
    dplyr::select(col_name, relevant_samples)
  colnames(fpkm_data)[-c(1:2)] <- sample_overview %>%
    dplyr::filter(sample %in% colnames(fpkm_data)[-c(1:2)]) %>%
    unite("sample_id", sample:tissue) %>%
    pull(sample_id)

  print(names(fpkm_data))
  print(names(deseq2_table))
  joined_df <- join_tables(deseq2_table, fpkm_data)


  joined_df <- joined_df %>% dplyr::mutate(overexpressed_in = ifelse(logFoldChange > 0,
    c_group[1],
    c_group[2]
  ))
  output_tb <- joined_df %>% dplyr::filter(padj < pval_threshold & abs(logFoldChange) > lfc_threshold)
  # Shorten names to less than 31 chars or excel gets angry

  output_tb
}

output_list <- purrr::map2(diffexp_tb_path, contrast_groups, build_output_table, fpkm = fpkm, col_name = col_name)
print(contrast_names)
names(output_list) <- purrr::map_chr(
  names(contrast_names),
  function(comp_name) {
    comp_name <- ifelse(nchar(comp_name) > 31,
      stringr::str_trunc(comp_name, 30, "right"),
      comp_name
    )
    comp_name
  }
)
print(names(output_list))
# names(output_list) <- names(contrast_names)
openxlsx::write.xlsx(x = output_list, file = output_path, )

library(DESeq2)
library(readr)
library(viper)
library(progeny)
library(RNAscripts)
library(magrittr)
library(decoupleR)

if (exists("snakemake")) {
  dds_obj_rlog <- snakemake@input[["dds_obj"]]
  fpkm_path <- snakemake@input[["fpkm"]]
  outpath <- snakemake@output[["sample_progeny_table"]]
  threads <- snakemake@threads
} else {
  BASE_PATH <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/cre_2022"
  dds_obj_rlog <- file.path(BASE_PATH, "deseq2/rlog_transform.RDS.gz")
  fpkm_path <- file.path(BASE_PATH, "fpkm/all.tsv")
  threads <- 1
}


## Read and import data
count_df_rlog <- readRDS(dds_obj_rlog) %>%
  assay() %>%
  as.matrix()

fpkm <- read_tsv(fpkm_path)

# Rename rownames
print(rownames(count_df_rlog))
count_df_rlog <- rename_count_rownames(
  count_table = count_df_rlog,
  fpkm_table = fpkm
)

# Import progeny network and run decoupler wmean
organism <- "mouse"

prog_net <- decoupleR::get_progeny(organism = organism, top = 100)

PathwayActivity_counts <- decoupleR::run_wmean(count_df_rlog,
  network = prog_net,
  .mor = "weight",
  times = 1000
) %>%
  dplyr::filter(statistic == "norm_wmean") %>%
  dplyr::mutate(prog_score = (1 - p_value) * sign(score)) %>%
  tidyr::pivot_wider(
    id_cols = "condition", names_from = "source",
    values_from = "prog_score"
  ) %>%
  tibble::column_to_rownames("condition") %>%
  as.matrix()


PathwayActivity_counts <- as.data.frame(t(PathwayActivity_counts))
PathwayActivity_counts$pathways <- rownames(PathwayActivity_counts)

readr::write_csv(PathwayActivity_counts, outpath)

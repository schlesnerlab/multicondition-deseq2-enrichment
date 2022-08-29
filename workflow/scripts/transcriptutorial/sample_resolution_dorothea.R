library(readr)
library(viper)
library(DESeq2)
library(dorothea)
library(RNAscripts)
library(magrittr)

# put whatever is your working directory here
if (exists("snakemake")) {
  dds_obj_rlog <- snakemake@input[["dds_obj"]]
  fpkm_path <- snakemake@input[["fpkm"]]
  outpath <- snakemake@output[["sample_dorothea_table"]]
  threads <- snakemake@threads
} else {
  BASE_PATH <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/apelin_exp"
  dds_obj_rlog <- file.path(BASE_PATH, "deseq2/rlog_transform.RDS.gz")
  fpkm_path <- file.path(BASE_PATH, "fpkm/all.tsv")
  threads <- 1
}



# count_df_rlog <- as.data.frame(read_csv("data/count_df_rlog.csv"))
count_df_rlog <- readRDS(dds_obj_rlog) %>%
  assay()
# row.names(count_df_rlog) <- count_df_rlog$gene

# count_df_rlog <- count_df_rlog[,-1]
# count_df_rlog <- count_df_rlog[complete.cases(count_df_rlog),]

fpkm <- read_tsv(fpkm_path)

rename_count_rownames <- function(count_table, fpkm_table) {
  filer <- fpkm_table %>%
    dplyr::filter(gene %in% rownames(count_table)) %>%
    dplyr::filter(!duplicated(gname))

  count_table <- count_table[filer$gene, ]
  rownames(count_table) <- filer$gname
  count_table
}

count_df_rlog <- rename_count_rownames(
  count_table = count_df_rlog,
  fpkm_table = fpkm
) 

#regulons <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B", "C"))
doro_net <-  decoupleR::get_dorothea(organism = "mouse",levels =  c("A", "B", "C"))

TF_activity <- decoupleR::run_wmean(count_df_rlog, network = doro_net, times = 1000, minsize = 5) %>% 
  dplyr::filter(statistic == "norm_wmean") %>% dplyr::select(source, condition, score) %>%
  tidyr::pivot_wider(names_from = condition, values_from = score)
                                    
#TF_activity <- run_viper(count_df_rlog,
#  regulons = regulons,
#  options = list(
#    method = "scale", minsize = 4,
#    eset.filter = FALSE, cores = threads,
#    verbose = FALSE
#  )
#) %>% as.data.frame()
### Preparing dorothea

## Dorothea/viper

#TF_activity$TF <- row.names(TF_activity)


readr::write_csv(TF_activity, outpath)

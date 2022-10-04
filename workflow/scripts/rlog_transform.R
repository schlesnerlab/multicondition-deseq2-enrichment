library(tidyverse)
library(DESeq2)
deseq_obj <- readRDS(snakemake@input[["dds_obj"]])

rld <- DESeq2::rlog(deseq_obj, blind = FALSE)

saveRDS(rld, snakemake@output[["rld"]])

norm_counts <- assay(rld, withDimnames = T) %>% as_tibble(rownames = "gene")
og_gene_names <- norm_counts$gene
norm_counts$gene <- stringr::str_extract(norm_counts$gene,
  pattern = "^ENS[A-Z0-9]*"
)
# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
species <- "mmusculus"
while (class(mart)[[1]] != "Mart") {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = glue::glue("{species}_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        stop()
      }
      # hop to next mirror
      mart <- switch(mart,
        useast = "uswest",
        uswest = "asia",
        asia = "www",
        www = {
          # wait before starting another round through the mirrors,
          # hoping that intermittent problems disappear
          Sys.sleep(30)
          "useast"
        }
      )
    }
  )
}


# df <- read.table(snakemake@input"]], sep='\t', header=1)
gene_name_type <- snakemake@config[["gene_name_type"]]
if (gene_name_type == "ENSEMBL") {
  g2g <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name"
    ),
    filters = "ensembl_gene_id",
    values = stringr::str_extract(norm_counts$gene,
      pattern = "^ENS[A-Z0-9]*"
    ),
    mart = mart,
  )
  colnames(g2g) <- c("gene", "gname")

  annotated <- dplyr::left_join(norm_counts, g2g)
  annotated$gname <- ifelse(annotated$gname == "" | is.na(annotated$gname),
    annotated$gene, annotated$gname
  )
} else if (gene_name_type == "ENTREZ_ID") {
  stop("to be implemented")
} else if (gene_name_type == "HGNC") {
  annotated <- norm_counts
  annotated$gname <- annotated$gene
}
# annotated$gene <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated$gene <- og_gene_names

readr::write_tsv(annotated, snakemake@output[["fpkm"]])

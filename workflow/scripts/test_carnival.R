####

library(progeny)

library(CARNIVAL)
library(OmnipathR)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(visNetwork)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(BiocParallel)
library(DESeq2)
library(decoupleR)
library(RNAscripts)
# library(biomaRt)
# Basic function to convert human to mouse gene names
# human = useMart(biomart= "ensembl", dataset = "hsapiens_gene_ensembl")
# mouse = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

mouse_human_homologs <- readr::read_tsv("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
  col_names = c("hgene", "hID", "mgene", "mID", "lcol")
)


if (exists("snakemake")) {
  dds_path <- snakemake@input[["dds_obj"]]
  diffexp_tb_path <- snakemake@input[["table"]]
  fpkm_path <- snakemake@input[["fpkm"]]
  carnival_output <- snakemake@output[["carnival_out"]]
  contrast_groups <- snakemake@wildcards[["contrast"]]
  s_groups <- snakemake@params[["s_groups"]]
  register(MulticoreParam(snakemake@threads))
  contrast_name <- contrast_groups
  plot_path <- snakemake@params[["plot_path"]]
  comp_groups <- snakemake@config[["signif_comp"]]
  color_scheme <- snakemake@config[["group_colors"]]
  cplex_path <- snakemake@config[["cplex_solver"]]
  stopifnot("Cplexpath doesnt exist please give path" = file.exists(cplex_path))
  temp_path <- snakemake@params[["temp_path"]]
  nr <- snakemake@resources[["nr"]]
  thread_num <- snakemake@threads
  mem_mb <- snakemake@resources[["mem_mb"]]
  time_limit <- (snakemake@resources[["time_min"]] - 20) * 60
  run_vanilla <- snakemake@params[["run_vanilla"]]
  perturbation_gene <- snakemake@params[["perturbation_gene"]]
  progeny_data <- "../data/progenyMembers.RData"
} else {
  BASE_ANALYSIS_DIR <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/APLNR_KO"
  dds_path <- file.path(paste0(BASE_ANALYSIS_DIR), "deseq2/all.rds")
  diffexp_tb_path <- file.path(
    paste0(BASE_ANALYSIS_DIR),
    "results/diffexp/condition/AplnKO_vs_basal.diffexp.tsv"
  )
  fpkm_path <- file.path(BASE_ANALYSIS_DIR, "fpkm/all.tsv")
  contrast_groups <- c("AplnKO", "basal")
  plot_path <- "."
  register(SerialParam())
  s_groups <- c("AplnKO", "basal")
  contrast_name <- glue::glue("{contrast_groups[[1]]} vs {contrast_groups[[2]]}")
  the_yaml <- yaml::read_yaml("../configs/VascAge_APLN_KO.yaml")
  comp_groups <- the_yaml$comp_groups
  color_scheme <- the_yaml$group_colors
  carnival_output <- "./test_output.RDS.gz"
  cplex_path <- "/home/heyer/software/external/CPLEX_Studio201/cplex/bin/x86-64_linux/cplex"
  thread_num <- 6
  temp_path <- "./"
  mem_mb <- 8192
  time_limit <- 3600
  run_vanilla <- TRUE
  progeny_data <- "./data/progenyMembers.RData"
}

# Read DESeq2 oobject and other tables
dds_obj <- readRDS(dds_path)
diffexp_tb <- read_tsv(diffexp_tb_path,
  col_names = c(
    "gene_id",
    "baseMean",
    "logFoldChange",
    "lfcSE", "stat",
    "pvalue", "padj"
  ),
  skip = 1
)
# Normalized_counts <- getVarianceStabilizedData(dds_obj)
Normalized_counts <- assay(vst(dds_obj, blind = F))
fpkm <- read_tsv(fpkm_path)


filer <- fpkm %>%
  dplyr::filter(gene %in% rownames(Normalized_counts)) %>%
  dplyr::filter(!duplicated(gname))

joined_df <- join_tables(diffexp_tb, filer) %>% dplyr::filter(!duplicated(gname))
Normalized_counts <- Normalized_counts[filer$gene, ]

rownames(Normalized_counts) <- filer$gname

joined_df %>%
  dplyr::select(gname, stat) %>%
  dplyr::filter(!is.na(stat)) %>%
  column_to_rownames(var = "gname") %>%
  as.matrix() -> diffexp_matrix




# regulons <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B"))
organism <- "mouse"
doro_net <- decoupleR::get_dorothea(organism = organism, levels = c("A", "B", "C"))
prog_net <- decoupleR::get_progeny(organism = organism, top = 100)

PathwayActivity <- PathwayActivity_CARNIVALinput <- run_wmean(
  mat = diffexp_matrix, network = prog_net, .source = "source", .target = "target",
  .mor = "weight", times = 1000, minsize = 5
) %>%
  dplyr::filter(statistic == "wmean") %>%
  as.data.frame()
if (any(abs(PathwayActivity$score) > 1)) {
  PathwayActivity$score <- sign(PathwayActivity$score) * (1 - PathwayActivity$p_value)
  warning("decoupler based enriched failed, falling back on (1-pvalue) * sign(score)")
}
# PathwayActivity <- PathwayActivity_CARNIVALinput <- progeny(diffexp_matrix,
#  scale = TRUE, organism = "Mouse", top = 100, perm = 10000, z_scores = F
# ) %>%
#  t() %>%
#  as.data.frame() %>%
#  tibble::rownames_to_column(var = "Pathway")

# colnames(PathwayActivity)[2] <- "score"
progeny_key <- setNames(
  object = PathwayActivity$score,
  nm = PathwayActivity$source
)

prog_net %>%
  mutate(progeny_activity = recode(source, !!!progeny_key)) %>%
  mutate(carnival_score = sign(weight) * progeny_activity) -> prog_net
prog_net %>%
  group_by(target) %>%
  summarise(carnival_score = mean(carnival_score)) -> prog_net_final

tf_activities_stat <- decoupleR::run_wmean(diffexp_matrix, network = doro_net, times = 1000, minsize = 5) %>%
  filter(statistic == "norm_wmean")
# options = list(
#    minsize = 5, eset.filter = FALSE,
#    cores = 1, verbose = FALSE, nes = TRUE
#  )
# )

prog_net_final %>% filter(!(target %in% tf_activities_stat$source)) -> prog_net_final

tf_activities <- tf_activities_CARNIVALinput <- tf_activities_stat %>% dplyr::select(source, score)


## Get Omnipath
omniR <- import_omnipath_interactions(organism = 10090)
# omniR <- import_pathwayextra_interactions(organism = 10090)
# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
  (consensus_stimulation == 1 |
    consensus_inhibition == 1
  ))

# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which(omnipath_sd$consensus_stimulation == 0)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 1)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 0)] <- 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[, c(
  "source_genesymbol", "consensus_stimulation",
  "consensus_inhibition", "target_genesymbol"
)] %>%
  dplyr::filter(consensus_stimulation == consensus_inhibition) %>%
  unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c("source", "interaction", "target")

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)

# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$source
tf_activities_carnival$source <- NULL
tf_list <- generateTFList(tf_activities_carnival, top = 50, access_idx = 1)
tf_vec <- tf_list$score[, ]

# For mouse change trp53 to tp53 and trp63 to tp63

names(tf_vec) <- gsub("^Trp", "Tp", names(tf_vec))

# progeny for CARNIVAL
# load(file = progeny_data)
# progenyMembers$gene$p53 <- "TP53"
# progmem_mouse <- purrr::map(progenyMembers$gene, convertHumanGeneHomologs, jax_database = mouse_human_homologs)
# progenyMembers$gene <- progmem_mouse
# progenyMembers$gene$p53 <- "Tp53"


# PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
# rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
# PathwayActivity_carnival$Pathway <- NULL
# progenylist <- assignPROGENyScores(
#  progeny = t(PathwayActivity_carnival),
#  progenyMembers = progenyMembers,
#  id = "gene",
#  access_idx = 1
# )
# progeny_vec <- progenylist$score
progeny_vec <- setNames(prog_net_final$carnival_score, nm = prog_net_final$target)
# get initial nodes

lp_opts <- CARNIVAL::defaultCplexCarnivalOptions(
  solverPath = cplex_path,
  cplexMemoryLimit = mem_mb,
  threads = thread_num,
  timelimit = time_limit,
  lpFilename = file.path(temp_path, "lptest.lp"),
  outputFolder = file.path(temp_path, "carnout")
)
cplex_opts <- suggestedCplexSpecificOptions()
lp_opts[names(cplex_opts)] <- cplex_opts
# lp_opts$solverPath <- cplex_path
lp_opts$cplexMemoryLimit <- mem_mb
lp_opts$threads <- thread_num
lp_opts$timelimit <- time_limit
# lp_opts$lpFilename <- file.path(temp_path, "lptest.lp")
# lp_opts$outputFolder <- file.path(temp_path, "carnout")
dir.create(lp_opts$outputFolder, showWarnings = F, recursive = T)

# run carnival
dir.create(file.path(temp_path, "carnout"), recursive = T, showWarnings = F)
# setwd(file.path(temp_path, "carnout"))
if (run_vanilla) {
  carnival_result <- runVanillaCarnival(
    perturbations = c(perturbation_gene = 1),
    measurements = unlist(tf_vec),
    priorKnowledgeNetwork = sif,
    weights = unlist(progeny_vec),
    carnivalOptions = lp_opts
  )
} else {
  carnival_result <- runInverseCarnival(
    measurements = unlist(tf_vec),
    priorKnowledgeNetwork = sif,
    weights = unlist(progeny_vec),
    carnivalOptions = lp_opts
  )
}
if (length(carnival_result$sifAll) < 50) {
  if (nr >= 2) {
    write(paste0("Failed to converge after 2 attempts"), file = stderr())
    carnival_result <- list()
  } else {
    stop(glue::glue("Attempt {nr}. Failed to converge restarting"))
  }
}
saveRDS(carnival_result, file = carnival_output)

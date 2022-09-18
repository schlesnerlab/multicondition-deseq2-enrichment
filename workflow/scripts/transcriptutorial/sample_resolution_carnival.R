library(readr)
library(CARNIVAL)
if (!require(RNAscripts)) {
  devtools::install("../scripts/RNAscripts", upgrade = "never", quiet = TRUE)
}
library(RNAscripts)
library(magrittr)
library(dplyr)
library(OmnipathR)
# files
if (exists("snakemake")) {
  tf_act_path <- snakemake@input[["tf_activities"]]
  path_act_path <- snakemake@input[["pathway_activities"]]
  cplex_path <- snakemake@config[["cplex_solver"]]
  carnival_sample_result <- snakemake@output[["carnival_sample_result"]]
  thread_num <- snakemake@threads
  time_limit <- (snakemake@resources[["time_min"]] - 20) * 60
  mem_mb <- snakemake@resources[["mem_mb"]]
  run_vanilla <- snakemake@params[["run_vanilla"]]
  temp_path <- snakemake@params[["temp_path"]]
  nr <- snakemake@resources[["nr"]]
  sample_id <- snakemake@wildcards[["sample"]]
  progeny_data <- "./workflow/data/progenyMembers.RData"
} else {
  base_path <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/APLN_KO"
  tf_act_path <- file.path(base_path, "dorothea/TF_act_sample_resolution.csv")
  path_act_path <- file.path(base_path, "progeny/sample_progney.csv")
  carnival_output <- "./test_output.RDS.gz"
  cplex_path <- "/home/heyer/software/external/CPLEX_Studio201/cplex/bin/x86-64_linux/cplex"
  thread_num <- 4
  temp_path <- file.path(base_path, "results/carnival/temp/VascAge_Young-EC-Plus-5")
  mem_mb <- 8192
  time_limit <- 3600
  run_vanilla <- FALSE
  sample_id <- "VascAge_Young-EC-Plus-5"
  progeny_data <- "./workflow/data/progenyMembers.RData"
  nr <- 1
}
if (nr == 4) {
  saveRDS(list(), file = carnival_sample_result)
  write("Failed to converge", file = stderr())
  quit(save = "no")
}
sample_id <- stringr::str_replace_all(sample_id, pattern = "-", replacement = ".")

tf_activities <- read_csv(tf_act_path)
pathway_activity <- read_csv(path_act_path)

mouse_human_homologs <- readr::read_tsv("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
  col_names = c("hgene", "hID", "mgene", "mID", "lcol")
)



## Get Omnipath
get_organism_number <- function(organism_name) {
  if (organism_name == "human") {
    organism_number = 9606
  } else if (organism_name == "mouse") {
    organism_number = 10090
  }
  organism_number
}
organism <- "mouse"
org_nb <- get_organism_number(organism)
omniR <- import_omnipath_interactions(organism = org_nb)


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
tfList <- generateTFList(tf_activities_carnival,
  top = 50,
  access_idx = 1:ncol(tf_activities_carnival)
)

# progeny for CARNIVAL
load(file = progeny_data)
#progenyMembers$gene$p53 <- "TP53"
if(organism == "mouse") {
  progmem_mouse <- purrr::map(progenyMembers$gene, convertHumanGeneHomologs, mouse_human_homologs)
  progenyMembers$gene <- progmem_mouse
  progenyMembers$gene$p53 <- "Tp53"
  names(progenyMembers$gene) [names(progenyMembers$gene) == "JAK.STAT"] <- "JAK-STAT"
} else {
  progenyMembers$gene$p53 <- "TP53"  
}
pathway_activity_carnival <- data.frame(pathway_activity, stringsAsFactors = F)
rownames(pathway_activity_carnival) <- pathway_activity_carnival %>% dplyr::pull(pathways)
# pathway_activity_carnival <- pathway_activity_carnival[,-c()]
pathway_activity_carnival$pathways <- NULL
progenylist <- assignPROGENyScores(
  progeny = t(pathway_activity_carnival),
  progenyMembers = progenyMembers,
  id = "gene",
  access_idx = 1:ncol(pathway_activity_carnival)
)
# get initial nodes


# run carnival for all samples
dir.create(file.path(temp_path, "carnout"), recursive = T, showWarnings = F)
# setwd(file.path(temp_path, "carnout"))
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
if (run_vanilla) {
  carnival_result <- runVanillaCarnival(
    perturbations = c("Aplnr" = 1, "Apln" = 1),
    measurements = unlist(tfList[[sample_id]]),
    priorKnowledgeNetwork = sif,
    weights = unlist(progenylist[[sample_id]]),
    carnivalOptions = lp_opts
  )
} else {
  carnival_result <- runInverseCarnival(
    measurements = unlist(tfList[[sample_id]]),
    priorKnowledgeNetwork = sif,
    weights = unlist(progenylist[[sample_id]]),
    carnivalOptions = lp_opts
  )

  # transoform to data.frame
  #  carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF,
  #  stringsAsFactors = F)
  #  carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
  #  carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)

  #  carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
  #  carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
  #  carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
  #  carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
  #  carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)
}
if (length(carnival_result$sifAll) < 50) {
  if (nr >= 2) {
    write(paste0("Failed to converge after 2 attempts"), file = stderr())
    carnival_result <- list()
  } else {
    stop("Failed to converge restarting")
  }
}

saveRDS(carnival_result, carnival_sample_result)

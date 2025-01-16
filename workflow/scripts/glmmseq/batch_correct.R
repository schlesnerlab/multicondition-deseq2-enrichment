library(sva)
library(readr)
library(DESeq2)

# snakemake Boilerplate
if (exists("snakemake", inherits = FALSE)) {
    counts <- snakemake@input[["counts"]]
    samples <- snakemake@params[["samples"]]
    batch_variable <- snakemake@params[["batch_variable"]]
    batch_correct_model <- snakemake@params[["batch_correct_model"]] %>% as.formula()
    batch_corrected_counts <- snakemake@output[["batch_corrected_counts"]]
    uncorrected_counts <- snakemake@output[["uncorrected_counts"]]
} else {
    counts <- "/desktop-home/heyer/projects/Vascular_Aging/Integrative_Analysis/join_RNA/count.tsv"
    samples <- "/desktop-home/heyer/projects/Vascular_Aging/Integrative_Analysis/join_RNA/full_metadata.tsv"
    batch_correct_model <- ~as.factor(age) + as.factor(Aplnr_KO) + as.factor(Apln_treatment) + as.factor(EC_status)
    batch_variable <- "experiment"  # Define a default batch variable if not using snakemake
    batch_corrected_counts <- "corrected_counts.rds"  # Define a default output path if not using snakemake
    uncorrected_counts <- "/desktop-home/heyer/projects/Vascular_Aging/Integrative_Analysis/join_RNA/dds_test.rds.gz"
}

# Read counts matrix
count_mat <-read.table(counts,  header = TRUE,
                       row.names = 1,
                       check.names = FALSE
)

sample_mat <- readr::read_tsv(samples)

# Ensure all count_mat columns are present in sample_mat
stopifnot(all(colnames(count_mat) %in% sample_mat$sample))
# Reorder sample mat to match count mat
sample_mat <- sample_mat[sample_mat$sample %in% colnames(count_mat),]

# Check if all variables in model are present in metadata
stopifnot(all(all.vars(batch_correct_model) %in% colnames(sample_mat)))

m <- model.matrix(batch_correct_model, data = sample_mat)

uncorrected_dds <- DESeqDataSetFromMatrix(count_mat, colData = sample_mat, design = ~condition)
write_rds(uncorrected_dds, uncorrected_counts )

corrected_counts <- sva::ComBat_seq(count_mat, batch = sample_mat |> dplyr::pull(!!batch_variable), covar_mod = m)

# Initialize DESeq2 object
corrected_dds <- DESeqDataSetFromMatrix(corrected_counts, colData = sample_mat, design = ~condition)
keep <- rowSums(counts(corrected_dds)) >= 30
corrected_dds <- corrected_dds[keep,]
corrected_dds <- DESeq(corrected_dds)

# Write corrected counts matrix
saveRDS(corrected_dds, file = batch_corrected_counts)
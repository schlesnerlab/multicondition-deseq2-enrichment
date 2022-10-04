### Snakefile
# Contains input functions and other functions for snakemake rules


def get_diffxp_files():
    """
    Return all files from DE and enrichment analyses.
    """
    output_files = []

    output_files.extend([join(BASE_ANALYSIS_DIR, "results/pca.svg")])
    output_files.extend(
        [
            join(BASE_ANALYSIS_DIR, "dorothea/dorothea_results.html"),
        ]
    )
    for cond in CONTRASTS:

        output_files.extend(
            expand(
                [
                    join(
                        BASE_ANALYSIS_DIR,
                        "reports/deseq2/{condition}_cohort_wide_comparison.html",
                    ),
                    join(
                        BASE_ANALYSIS_DIR,
                        "reports/deseq2/{condition}_joint_gsea_report.html",
                    ),
                    join(
                        BASE_ANALYSIS_DIR,
                        "results",
                        "diffexp",
                        "excel_tables",
                        "{condition}_diffexp_genes.xlsx",
                    ),
                ],
                condition=cond,
            )
        )
        output_files.extend(
            expand(
                [
                    join(
                        BASE_ANALYSIS_DIR,
                        "results/diffexp/{condition}/{contrast}.diffexp.tsv",
                    ),
                    join(
                        BASE_ANALYSIS_DIR,
                        "results/diffexp/{condition}/{contrast}.ma-plot.pdf",
                    ),
                    join(
                        BASE_ANALYSIS_DIR,
                        "reports/deseq2/{condition}/{contrast}_diffexp.html",
                    ),
                    join(
                        BASE_ANALYSIS_DIR, "progeny/{condition}/{contrast}_progeny.html"
                    ),
                    join(
                        BASE_ANALYSIS_DIR,
                        "dorothea/{condition}/{contrast}_dorothea.html",
                    ),
                ],
                contrast=config["diffexp"]["contrasts"][cond],
                condition=cond,
            )
        )
    return output_files


def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


### diffexp.smk
def get_count_matrix(wildcards):
    """
        Gets the path to the count matrix file. depending on the Infrastructure used.
    TODO: Expand once we have new functions for expression data in place
    """
    if "counts" in config:
        count_file = config["counts"]
    else:
        count_file = join(BASE_ANALYIS_DIR, "counts/all.tsv")
    return count_file


def get_deseq2_threads(wildcards=None):
    """
    https://twitter.com/mikelove/status/918770188568363008
    """
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.condition][wildcards.contrast]


def get_all_contrasts(wildcards):
    return config["diffexp"]["contrasts"][wildcards.condition]


def get_gsea_results(wildcards):
    """
    Helper function to retrieve fitting gsea results
    """
    cond = wildcards.condition
    gsea_files = expand(
        join(BASE_ANALYSIS_DIR, "results/diffexp/{condition}/{contrast}.gseares.RDS"),
        contrast=CONTRASTS[cond],
        condition=cond,
    )
    return gsea_files


def get_diffexp_tables(wildcards):
    """
    Retrieves all DE results for one condition given.
    """
    cond = wildcards.condition
    output_files = expand(
        join(BASE_ANALYSIS_DIR, "results/diffexp/{{condition}}/{contrast}.diffexp.tsv"),
        contrast=config["diffexp"]["contrasts"][cond],
    )
    return output_files


def get_carnival_objs(wildcards):
    """
    Retrieves Carnival results for one condition
    """
    cond = wildcards.condition
    output_files = expand(
        join(
            BASE_ANALYSIS_DIR,
            "results/{{type}}/{{condition}}/{contrast}_carnival_res.RDS.gz",
        ),
        contrast=config["diffexp"]["contrasts"][cond],
    )
    return output_files


def get_mem_mb(wildcards, attempt):
    return (attempt * 20480) + 40960


def get_sample(wildcards):
    return samples.loc[wildcards.sample]["sample"]

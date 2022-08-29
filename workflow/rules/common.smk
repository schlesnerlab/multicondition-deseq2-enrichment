### Snakefile


def get_diffxp_files():
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
                        "results/diffexp/{condition}/{contrast}.ma-plot.svg",
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
    Gets the path to the count matrix file. depending on the INfrastructure used.
    """
    if "counts" in config:
        return config["counts"]
    else:
        return join(BASE_ANALYIS_DIR, "counts/all.tsv")


def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.condition][wildcards.contrast]


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
    cond = wildcards.condition
    output_files = expand(
        join(BASE_ANALYSIS_DIR, "results/diffexp/{{condition}}/{contrast}.diffexp.tsv"),
        contrast=config["diffexp"]["contrasts"][cond],
    )
    return output_files


### CARNIVAL
def get_mem_mb(wildcards, attempt):
    return (attempt * 20480) + 40960


def get_sample(wildcards):
    return samples.loc[wildcards.sample]["sample"]

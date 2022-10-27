rule count_matrix:
    input:
        expand(
            join(
                BASE_DATA_DIR,
                "{unit.sample}",
                "{unit.tissue}",
                SEQ_TYPE,
                "merged-alignment/featureCounts/{unit.tissue}_{unit.sample}.fpkm_tpm.featureCounts.tsv",
            ),
            unit=samples.itertuples(),
        ),
    output:
        join(BASE_ANALYSIS_DIR, "counts/all.tsv"),
    params:
        samples=samples["sample"].tolist(),
        strand=get_strandness(samples),
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/count_matrix.log",
    script:
        "../scripts/count-matrix.py"

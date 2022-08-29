

rule run_carnival_vanilla:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        table=join(BASE_ANALYSIS_DIR, "results/diffexp/{contrast}.diffexp.tsv"),
    output:
        carnival_out=join(
            BASE_ANALYSIS_DIR, "results/carnival/{contrast}_carnival_res.RDS.gz"
        ),
    params:
        s_groups=samples.condition.unique(),
        temp_path=temp(join(BASE_ANALYSIS_DIR, "results/carnival/temp/{contrast}/")),
        run_vanilla=True,
    conda:
        "../envs/carnival.yaml"
    threads: 12
    resources:
        nr=lambda wildcards, attempt: attempt,
        mem_mb=81920,
        time_min=(9 * 60) + 20,
    log:
        "logs/carnival/{contrast}_carn.log",
    script:
        "../scripts/test_carnival.R"


rule run_inverse_carnival:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/all.rds"),
        table=join(BASE_ANALYSIS_DIR, "results/diffexp/{contrast}.diffexp.tsv"),
    output:
        carnival_out=join(
            BASE_ANALYSIS_DIR, "results/inversecarnival/{contrast}_carnival_res.RDS.gz"
        ),
    params:
        s_groups=samples.condition.unique(),
        temp_path=temp(
            join(BASE_ANALYSIS_DIR, "results/inversecarnival/temp/{contrast}/")
        ),
        run_vanilla=False,
    conda:
        "../envs/carnival.yaml"
    threads: 20
    resources:
        mem_mb=get_mem_mb,
        time_min=(25 * 60) + 40,
        nr=lambda wildcards, attempt: attempt,
    log:
        "logs/carnival/{contrast}_carn.log",
    script:
        "../scripts/test_carnival.R"


rule carnival_deseq_report:
    input:
        carnival_obj=join(
            BASE_ANALYSIS_DIR, "results/{type}/{contrast}_carnival_res.RDS.gz"
        ),
    output:
        carnival_rep=join(
            BASE_ANALYSIS_DIR, "results/{type}/reports/{contrast}_results.html"
        ),
    params:
        tutorial_source_path="scripts/transcriptutorial/",
    conda:
        "../envs/carnival.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=59,
    log:
        "logs/carnival/{type}_{contrast}_result.log",
    script:
        "../scripts/RMD_scripts/carnival_downstream.Rmd"


rule carnival_joint_report:
    input:
        carnival_objs=expand(
            join(BASE_ANALYSIS_DIR, "results/{{type}}/{contrast}_carnival_res.RDS.gz"),
            contrast=config["diffexp"]["contrasts"],
        ),
    output:
        carnival_rep=join(BASE_ANALYSIS_DIR, "results/{type}/reports/join_report.html"),
    params:
        tutorial_source_path="scripts/transcriptutorial/",
    conda:
        "../envs/carnival.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=20,
    log:
        "logs/carnival/{type}_report.log",
    script:
        "../scripts/RMD_scripts/carnival_join.Rmd"


rule run_sample_dorothea:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/rlog_transform.RDS.gz"),
    output:
        sample_dorothea_table=join(
            BASE_ANALYSIS_DIR, "dorothea/TF_act_sample_resolution.csv"
        ),
    conda:
        "../envs/carnival.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=59,
    log:
        "logs/sample_dorothea/report.log",
    script:
        "../scripts/transcriptutorial/sample_resolution_dorothea.R"


rule run_sample_progeny:
    input:
        dds_obj=join(BASE_ANALYSIS_DIR, "deseq2/rlog_transform.RDS.gz"),
    output:
        sample_progeny_table=join(BASE_ANALYSIS_DIR, "progeny/sample_progney.csv"),
    conda:
        "../envs/carnival.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=59,
    log:
        "logs/sample_dorothea/report.log",
    script:
        "../scripts/transcriptutorial/sample_resolution_progeny.R"


rule run_sample_carnival:
    input:
        tf_activities=join(BASE_ANALYSIS_DIR, "dorothea/TF_act_sample_resolution.csv"),
        pathway_activities=join(BASE_ANALYSIS_DIR, "progeny/sample_progney.csv"),
    output:
        carnival_sample_result=join(
            BASE_ANALYSIS_DIR,
            "results/carnival/sample_resolution/{sample}_carn_res.RDS.gz",
        ),
    params:
        run_vanilla=False,
        sample_id=get_sample,
        temp_path=temp(join(BASE_ANALYSIS_DIR, "results/carnival/temp/{sample}/")),
    conda:
        "../envs/carnival.yaml"
    threads: 8
    resources:
        mem_mb=get_mem_mb,
        time_min=(9 * 60) + 20,
        nr=lambda wildcards, attempt: attempt,
    log:
        "logs/sample_dorothea/{sample}_report.log",
    script:
        "../scripts/transcriptutorial/sample_resolution_carnival.R"


rule sample_carnival_report:
    input:
        carnival_sample_result=expand(
            join(
                BASE_ANALYSIS_DIR,
                "results/carnival/sample_resolution/{sample.sample}_carn_res.RDS.gz",
            ),
            sample=samples.itertuples(),
        ),
    output:
        join(BASE_ANALYSIS_DIR, "reports/carnival/carnival_sample_report.html"),
    params:
        tutorial_source_path="scripts/transcriptutorial/",
        sample_names=samples["sample"],
    conda:
        "../envs/carnival.yaml"
    threads: 1
    resources:
        mem_mb=8192,
        time_min=120,
    log:
        "logs/sample_dorothea/report.log",
    script:
        "../scripts/transcriptutorial/06_analysis_CARNIVAL_results.Rmd"

rule batch_correct:
    input: 
        counts = get_count_matrix,
    output: 
        batch_corrected_counts =  join(BASE_ANALYSIS_DIR, "counts/batch_corrected_counts.rds")
    params:
        samples = config["samples"],
        batch_correct_model = config["batch_correct"]["model"]
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/batch_correct.log"
    resources:
        mem_mb=8192,
        time_min=59
    threads: 4
    script:
        "../scripts/glmmseq/batch_correct.R"

rule run_glmmseq:
    input:
        batch_corrected_counts = join(BASE_ANALYSIS_DIR, "counts/batch_corrected_counts.rds")
    output:
        glmmseq_obj = join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_obj.rds.gz"),
        glmmseq_refit = join(BASE_ANALYSIS_DIR, "glmmseq/glmmseq_refit.rds.gz"),
    params:
        formula = config["glmmseq"]["formula"],
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/glmmseq/run_glmmseq.log"
    resources:
        mem_mb = 16384,
        time_min = 59
    threads: 12
    script:
        "../scripts/glmmseq/glmmseq.R"


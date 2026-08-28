wildcard_constraints:
    metric = "|".join(JUI_METRICS.keys()),

rule ComputeJUI:
    """
    Per-sample Junction Usage Index (my_utils' compute-jui CLI). Purely per-sample -- no
    expand() here -- so adding a sample never retriggers existing samples' JUI computation,
    unlike leafcutter's clustering step. --decompose is always on so every downstream matrix
    flavor (usage/IR/AS) in config/jui_metrics.yaml is derivable from this one file.
    """
    input:
        junc = "SplicingAnalysis/juncfiles/{sample}.junc",
        bigwig = "bigwigs/unstranded_raw/{sample}.bw",
    params:
        plain_tsv = "SplicingAnalysis/JUI/{sample}.jui.tsv",
    output:
        tsv = "SplicingAnalysis/JUI/{sample}.jui.tsv.gz",
        tbi = "SplicingAnalysis/JUI/{sample}.jui.tsv.gz.tbi",
    conda:
        "../envs/jui.yml"
    log:
        "logs/ComputeJUI/{sample}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(8000, 16000)
    shell:
        """
        compute-jui --junc-table {input.junc} --junc-format bed12 \
            --bigwig {input.bigwig} --strand unstranded --decompose \
            --output-tsv {params.plain_tsv} --tabix &> {log}
        """

rule JUI_MergeMatrix:
    """
    Pool per-sample JUI tsvs into one junctions x samples matrix, per GenomeName (samples on
    different genomes are never pooled together, matching leafcutter_cluster's convention).
    Unlike ComputeJUI, this rule genuinely reruns in full whenever a sample is added -- the
    merge engine (scripts/MergeJUIMatrix.py) is a streaming heapq k-way merge chosen specifically
    to keep that rerun's memory flat regardless of cohort size (see the project plan for the
    benchmark that motivated this over DuckDB/pandas/polars).
    """
    input:
        tsvs = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("SplicingAnalysis/JUI/{sample}.jui.tsv.gz"),
    params:
        expr = lambda wildcards: JUI_METRICS[wildcards.metric]["expr"],
        fill = lambda wildcards: JUI_METRICS[wildcards.metric]["fill"],
    output:
        matrix = "SplicingAnalysis/JUI/{GenomeName}/matrix.{metric}.tsv.gz",
        tbi = "SplicingAnalysis/JUI/{GenomeName}/matrix.{metric}.tsv.gz.tbi",
    conda:
        "../envs/jui.yml"
    log:
        "logs/JUI_MergeMatrix/{GenomeName}.{metric}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(2000, 4000)
    shell:
        """
        python scripts/MergeJUIMatrix.py --input-tsvs {input.tsvs} \
            --metric-expr {params.expr:q} --fill {params.fill} \
            --output {output.matrix} --tabix &> {log}
        """


rule featurecounts:
    input:
        bam = ExpandAllSTARSamplesInFormatStringFromGenomeNameAndStrandWildcards("Alignments/{sample}/Aligned.sortedByCoord.out.bam"),
        index = ExpandAllSTARSamplesInFormatStringFromGenomeNameAndStrandWildcards("Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
    output:
        counts = "featureCounts/{GenomeName}/{Strandedness}.Counts.txt",
        summary = "featureCounts/{GenomeName}/{Strandedness}.Counts.txt.summary",
    threads:
        8
    conda:
        "../envs/subread_featureCounts.yml"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 48000),
        tasks = 9,
    log:
        "logs/featureCounts/{GenomeName}.{Strandedness}.log"
    params:
        strand = lambda wildcards: {'FR':'-s 1', 'U':'-s 0', 'RF':'-s 2'}[wildcards.Strandedness],
        paired_flag = UsePairedEndFeatureCountsIfMixingSingleAndPairedReads,
        extra = ""
    shell:
        """
        featureCounts {params.strand} {params.paired_flag} {params.extra} -T {threads} --ignoreDup --primary -a {input.gtf} -o {output.counts} {input.bam} --maxMOp 100 &> {log}
        """

use rule featurecounts as featurecounts_allUnstranded with:
    input:
        bam = ExpandAllSTARSamplesInFormatStringFromGenomeNameWildcard("Alignments/{sample}/Aligned.sortedByCoord.out.bam"),
        index = ExpandAllSTARSamplesInFormatStringFromGenomeNameWildcard("Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
    output:
        counts = "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt",
        summary = "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt.summary",
    threads:
        8
    resources:
        mem_mb = 12000,
        tasks = 9,
    log:
        "logs/featureCounts/{GenomeName}.AllUnstranded.log"
    params:
        strand = '-s 0',
        paired_flag = UsePairedEndFeatureCountsIfMixingSingleAndPairedReads,
        extra = ""

rule edgeR_differential_expression:
    input:
        groupfile = lambda wildcards: os.path.abspath(config['contrast_group_files_prefix'] + "{contrast}.txt"),
        counts = lambda wildcards: get_filled_path_from_contrast(wildcards, "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt"),
    output:
        results = "differential_expression/{contrast}/results.tsv.gz",
        plots = directory("differential_expression/{contrast}/plots")
    conda:
        "../envs/r_2.yml"
    log:
        "logs/differential_expression/{contrast}.log"
    shell:
        """
        Rscript scripts/BasicDifferentialExpression_edgeR.R {input.counts} {input.groupfile} {output.results} {output.plots} &> {log}
        """

rule ExpressionMatrix:
    input:
        counts = "featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt",
    output:
        tpm = "ExpressionMatrices/{GenomeName}/log2TPM.bed",
        cpm = "ExpressionMatrices/{GenomeName}/log2TMM_CPM.bed",
        filtered_cpm = "ExpressionMatrices/{GenomeName}/log2Filtered_TMM_CPM.bed",
    conda:
        "../envs/rnanorm.yaml"
    log:
        "logs/ExpressionMatrix/{GenomeName}.log"
    params:
        sample_rename_regex = "'.+\/(.+?)\/Aligned\.sortedByCoord\.out\.bam'",
        extra = "--pseudocount 0.1"
    resources:
        mem_mb = 16000,
        tasks = 1,
    shell:
        """
        python scripts/CountsToExpressionMatrix.py -i {input.counts} --log2TPM_Matrix_bed {output.tpm} --log2TMM_Normalized_CPM_Matrix_bed {output.cpm} --log2Filtered_TMM_Normalized_CPM_Matrix_bed {output.filtered_cpm} --sample_rename_regex {params.sample_rename_regex} {params.extra} &> {log}
        """
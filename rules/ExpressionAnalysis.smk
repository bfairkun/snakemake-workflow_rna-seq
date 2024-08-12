
rule featurecounts:
    input:
        bam = ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam"),
        index = ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
    output:
        counts = "featureCounts/{GenomeName}/{Strandedness}.Counts.txt",
        summary = "featureCounts/{GenomeName}/{Strandedness}.Counts.txt.summary",
    threads:
        8
    conda:
        "../envs/subread_featureCounts.yml"
    resources:
        mem_mb = 12000,
        tasks = 9,
    log:
        "logs/featureCounts/{GenomeName}.{Strandedness}.log"
    params:
        strand = lambda wildcards: {'FR':'-s 1', 'U':'-s 0', 'RF':'-s 2'}[wildcards.Strandedness],
        extra = UsePairedEndFeatureCountsIfMixingSingleAndPairedReads
    shell:
        """
        featureCounts {params.strand} {params.extra} -T {threads} --ignoreDup --primary -a {input.gtf} -o {output.counts} {input.bam} &> {log}
        """

use rule featurecounts as featurecounts_allUnstranded with:
    input:
        bam = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam"),
        index = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
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
        extra = UsePairedEndFeatureCountsIfMixingSingleAndPairedReads

# rule GetGeneNames_bioMart:

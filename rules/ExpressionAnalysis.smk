
rule featurecounts:
    input:
        bam = ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam"),
        bai = ExpandAllSamplesInFormatStringFromGenomeNameAndStrandWildcards("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai"),
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
    output:
        "featureCounts/{GenomeName}/{Strandedness}.Counts.txt"
    threads:
        8
    resources:
        mem_mb = 12000,
        tasks = 9,
        # cpus_per_node = 9,
    log:
        "logs/featureCounts/{GenomeName}.{Strandedness}.log"
    params:
        strand = lambda wildcards: {'FR':'-s 1', 'U':'-s 0', 'RF':'-s 2'}[wildcards.Strandedness],
        extra = ""
    shell:
        """
        featureCounts {params.strand} {params.extra} -T {threads} --ignoreDup --primary -a {input.gtf} -o {output} {input.bam} &> {log}
        """


rule featurecounts:
    input:
        bam = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam"),
        bai = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai"),
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
    output:
        "featureCounts/{GenomeName}/Counts.txt"
    threads:
        8
    resources:
        mem = 12000,
        cpus_per_node = 9,
    log:
        "logs/featureCounts/{GenomeName}.log"
    params:
        extra = "-s 2"
    shell:
        """
        featureCounts {params.extra} -T {threads} --ignoreDup --primary -a {input.gtf} -o {output} {input.bam} &> {log}
        """

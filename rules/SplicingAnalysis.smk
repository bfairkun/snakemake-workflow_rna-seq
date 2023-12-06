rule ExtractJuncs:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        "SplicingAnalysis/juncfiles/{sample}.junc",
    params:
        strand = "0"
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{sample}.log"
    shell:
        """
        (regtools junctions extract -m 20 -s {params.strand} {input.bam} > {output}) &> {log}
        """


rule annotate_juncfiles:
    input:
        fa = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa"),
        fai = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai"),
        gtf = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"),
        juncs = "SplicingAnalysis/juncfiles/{sample}.junc",
    output:
        counts = "SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz"
    log:
        "logs/annotate_juncfiles/{sample}.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions annotate {input.juncs} {input.fa} {input.gtf} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{$4=$1"_"$2"_"$3"_"$6; print $4, $5}}' | gzip - > {output.counts} ) &> {log}
        """


rule ConcatJuncFilesAndKeepUniq:
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard("SplicingAnalysis/juncfiles/{sample}.junc")
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc"
    log:
        "logs/ConcatJuncFilesAndKeepUniq/{GenomeName}.log"
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/Collapse_Juncsfiles.R {output} {input} &> {log}
        """

rule AnnotateConcatedUniqJuncFile_basic:
    input:
        junc = "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa"
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    log:
        "logs/AnnotateConcatedUniqJuncFile_hg38Basic.{GenomeName}.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions annotate {input.junc} {input.fa} {input.gtf} | gzip - > {output} ) &> {log}
        """

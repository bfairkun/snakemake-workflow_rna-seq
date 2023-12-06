# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
# Config file should have unique sample names for each sample. Rows with the
# same samples name will be merged (eg, for combining the same sample
# sequencing across multiple lanes with multiple sets of fastq)
configfile: "config/config.yaml"

include: "rules/common.smk"

wildcard_constraints:
    GenomeName = "|".join(STAR_genomes.index),
    sample = "|".join(samples.index)
localrules: DownloadFastaAndGtf, CopyFastq, MultiQC

include: "rules/PreprocessAndAlign.smk"
include: "rules/IndexGenome.smk"
include: "rules/SplicingAnalysis.smk"
include: "rules/ExpressionAnalysis.smk"
include: "rules/QC.smk"
include: "rules/MakeBigwigs.smk"

rule all:
    input:
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",sample=samples.index),
        expand("SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz", sample=samples.index),
        "../output/QC/ReadCountsPerSamples.tsv",
        expand("bigwigs/unstranded/{sample}.bw", sample=samples.index),
        "SplicingAnalysis/ObservedJuncsAnnotations/GRCh38_GencodeRelease44Comprehensive.uniq.annotated.tsv.gz",
        "Multiqc/multiqc_report.html",
        "featureCounts/GRCh38_GencodeRelease44Comprehensive/Counts.txt"


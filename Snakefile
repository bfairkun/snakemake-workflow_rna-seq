# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
# Config file should have unique sample names for each sample. Rows with the
# same samples name will be merged (eg, for combining the same sample
# sequencing across multiple lanes with multiple sets of fastq)
configfile: "config/config.yaml"

include: "rules/common.py"

wildcard_constraints:
    GenomeName = "|".join(STAR_genomes.index),
    sample = "|".join(AllSamples),
    Strandedness = "|".join(["U", "FR", "RF"])

include: "rules/PreprocessAndAlign.smk"
include: "rules/IndexGenome.smk"
include: "rules/ExpressionAnalysis.smk"
include: "rules/SplicingAnalysis.smk" #Contains rules inherited from ExpressionAnalysis.smk
include: "rules/QC.smk"
include: "rules/MakeBigwigs.smk"

rule all:
    input:
        "samples.SRA_accession_links_filled.tsv",
        expand("Alignments/{sample}/Aligned.sortedByCoord.out.bam",sample=AllSamples),
        expand("SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz", sample=AllSamples),
        "../output/QC/ReadCountsPerSamples.tsv",
        expand("bigwigs/unstranded/{sample}.bw", sample=AllSamples),
        # "SplicingAnalysis/ObservedJuncsAnnotations/GRCh38_GencodeRelease44Comprehensive.uniq.annotated.tsv.gz",
        # "Multiqc",
        expand(config['GenomesPrefix'] + "{GenomeName}/STARIndex", GenomeName = samples['STARGenomeName'].unique()),
        expand("featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt", GenomeName = samples['STARGenomeName'].unique()),
        expand("SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.with_ss_scores.tsv.gz", GenomeName = samples['STARGenomeName'].unique()),
        expand("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/{Metric}.sorted.bed.gz", GenomeName = samples['STARGenomeName'].unique(), Metric = ["JuncCounts", "PSI", "PSI_ByMax", "PSIDenom", "PSI_ByMaxDenom"]),
        expand("SplicingAnalysis/ClassifyJuncs/{GenomeName}/Leaf2_junction_classifications.txt", GenomeName = samples['STARGenomeName'].unique()),
        expand("SplicingAnalysis/differential_splicing/{contrast}/", contrast=contrasts),
        expand("differential_expression/{contrast}/results.tsv.gz", contrast=contrasts),
        "Multiqc",
        # expand("QC/QualimapRnaseq/{sample}", sample=AllSamples),
        # config['GenomesPrefix'] + "GRCh38_GencodeRelease44Comprehensive/Reference.Transcripts.colored.bed.gz",
        # expand("featureCounts/GRCh38_GencodeRelease44Comprehensive/{Strandedness}.Counts.txt", Strandedness=samples['Strandedness'].unique())

rule Gather_Fastp_Fastqs:
    input:
        expand("FastqFastp/{sample}.fastp.html", sample=AllSamples),


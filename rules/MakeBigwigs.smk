
rule MakeBigwigs_Raw:
    """
    Unnormalized (raw) coverage bigwig. Backs JUI's intron-retention counts (rules/JUI.smk) and
    is the base for MakeBigwigs_NormalizedToGenomewideCoverage below, which rescales this rather
    than re-running samtools+bedtools genomecov over the whole BAM a second time.
    """
    input:
        fai = lambda wildcards: config['GenomesPrefix'] + samples.loc[samples['sample']==wildcards.sample]['STARGenomeName'].tolist()[0] + "/Reference.fa.fai",
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
    params:
        GenomeCovArgs="-split",
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
    shadow: "shallow"
    output:
        bw = "bigwigs/unstranded_raw/{sample}.bw",
        bw_minus = []
    log:
        "logs/MakeBigwigs_Raw/{sample}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(42000, 52000)
    conda:
        "../envs/pybedtools.yml"
    shell:
        """
        scripts/BamToBigwig.sh {input.fai} {input.bam} {output.bw}  GENOMECOV_ARGS="{params.GenomeCovArgs}" REGION='{params.Region}' MKTEMP_ARGS="{params.MKTEMP_ARGS}" SORT_ARGS="{params.SORT_ARGS}" {params.bw_minus}"{output.bw_minus}" &> {log}
        """

rule MakeBigwigs_NormalizedToGenomewideCoverage:
    """
    Scale the raw bigwig (MakeBigwigs_Raw) to coverage per billion covered bases genome-wide,
    reading/writing via pyBigWig (scripts/NormalizeBigwig.py) instead of a bam-derived read
    count -- total bases covered is the more principled denominator for scaling a coverage
    track, and reading it straight from the bigwig avoids a separate samtools idxstats pass.
    """
    input:
        bw_raw = "bigwigs/unstranded_raw/{sample}.bw",
    output:
        bw = "bigwigs/unstranded/{sample}.bw",
    log:
        "logs/MakeBigwigs_unstranded/{sample}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(4000, 8000)
    conda:
        "../envs/pybedtools.yml"
    shell:
        """
        python scripts/NormalizeBigwig.py --input {input.bw_raw} --output {output.bw} &> {log}
        """


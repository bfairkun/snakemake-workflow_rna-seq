
rule MakeBigwigs_NormalizedToGenomewideCoverage:
    """
    Scale bigwig to base coverage per billion chromosomal reads
    """
    input:
        fai = lambda wildcards: config['GenomesPrefix'] + samples.loc[samples['sample']==wildcards.sample]['STARGenomeName'].tolist()[0] + "/Reference.fa.fai",
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
        NormFactorsFile = "../output/QC/ReadCountsPerSamples.tsv"
    params:
        GenomeCovArgs="-split",
        bw_minus = "bw_minus=",
        MKTEMP_ARGS = "-p " + config['scratch'],
        SORT_ARGS="-T " + config['scratch'],
        Region = "",
    shadow: "shallow"
    output:
        bw = "bigwigs/unstranded/{sample}.bw",
        bw_minus = []
    log:
        "logs/MakeBigwigs_unstranded/{sample}.log"
    resources:
        mem_mb = much_more_mem_after_first_attempt
    shell:
        """
        ScaleFactor=$(bc <<< "scale=3;1000000000/$(grep '{wildcards.sample}' {input.NormFactorsFile} | awk 'NR==1 {{print $2}}')")
        scripts/BamToBigwig.sh {input.fai} {input.bam} {output.bw}  GENOMECOV_ARGS="{params.GenomeCovArgs} -scale ${{ScaleFactor}}" REGION='{params.Region}' MKTEMP_ARGS="{params.MKTEMP_ARGS}" SORT_ARGS="{params.SORT_ARGS}" {params.bw_minus}"{output.bw_minus}" &> {log}
        """


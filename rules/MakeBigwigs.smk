
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
        mem_mb = GetMemForSuccessiveAttempts(42000, 52000)
    conda:
        "../envs/pybedtools.yml"
    shell:
        """
        ReadCount=$(awk -F'\\t' -v f="idxstats/{wildcards.sample}.idxstats.txt" '$1==f {{print $2}}' {input.NormFactorsFile})
        if [ -z "$ReadCount" ]
        then
            echo "ERROR: no row for idxstats/{wildcards.sample}.idxstats.txt in {input.NormFactorsFile}" >&2
            exit 1
        fi
        ScaleFactor=$(bc <<< "scale=3;1000000000/$ReadCount")
        scripts/BamToBigwig.sh {input.fai} {input.bam} {output.bw}  GENOMECOV_ARGS="{params.GenomeCovArgs} -scale ${{ScaleFactor}}" REGION='{params.Region}' MKTEMP_ARGS="{params.MKTEMP_ARGS}" SORT_ARGS="{params.SORT_ARGS}" {params.bw_minus}"{output.bw_minus}" &> {log}
        """


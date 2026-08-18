rule QualimapRnaseq:
    input:
        gtf = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"),
        bam="Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        index="Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"
    output:
        results = "QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt",
        outdir = directory("QC/QualimapRnaseq/{sample}")
    log:
        "logs/QualimapRnaseq/{sample}.log"
    conda:
        "../envs/qualimap.yml"
    params:
        # extra = "-p strand-specific-reverse"
        extra = ""
    resources:
        mem_mb = 16000
    shell:
        """
        unset DISPLAY
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} {params.extra} --java-mem-size=12G -outdir {output.outdir}/ &> {log}
        """

rule MultiQC:
    input:
        expand("Alignments/{sample}/Log.final.out", sample= samples_STAR_aligned_here),
        expand("QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt", sample=AllSamples),
        expand("FastqFastp/{sample}.fastp.json", sample=samples_needing_alignment),
        expand("featureCounts/{GenomeName}/AllSamplesUnstrandedCounting.Counts.txt.summary", GenomeName=samples['STARGenomeName'].unique(), Strandedness=samples['Strandedness'].unique()),
        expand("idxstats/{sample}.idxstats.txt", sample=AllSamples)
    log: "logs/Multiqc.log"
    output:
        directory("Multiqc")
    localrule: True
    shell:
        """
        multiqc -f -o {output}/ {input} &> {log}
        """



rule CountReadsPerSample:
    input:
        expand("idxstats/{sample}.idxstats.txt", sample=AllSamples)
    output:
        "../output/QC/ReadCountsPerSamples.tsv"
    log:
        "logs/CountReadsPerSample.log"
    shell:
        """
        rm -f {output}.tmp
        for f in {input}
        do
           ReadCount=$(awk -F'\\t' '$1~"^chr[1-9]" {{sum+=$3}} END {{print sum+0}}' $f)
           if [ "$ReadCount" -eq 0 ]
           then
               echo "ERROR: $f has no reads on chromosomes matching ^chr[1-9], so bigwig normalization would be undefined. Check that this genome uses UCSC-style chromosome names." >&2
               exit 1
           fi
           printf "%s\\t%s\\n" $f $ReadCount >> {output}.tmp
        done
        mv {output}.tmp {output}
        """

rule CountMappedBasesPerSample:
    input:
        bam="Alignments/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "QC/MappedBases/{sample}.mappedBases.txt"
    log:
        "logs/CountMappedBasesPerSample.{sample}.log"
    shell:
        """
        bedtools bamtobed -split -bed12 -i {input.bam} | awk '{{
            n_blocks = $10
            split($11, block_sizes, ",")
            for (i = 1; i <= n_blocks; i++) {{
                if (block_sizes[i] != "") sum += block_sizes[i]
            }}
        }} END {{print "{wildcards.sample}\\t" sum}}' > {output}
        """

rule CollectMappedBasesPerSample:
    input:
        expand("QC/MappedBases/{sample}.mappedBases.txt", sample=AllSamples)
    output:
        "../output/QC/MappedBasesPerSamples.tsv"
    log:
        "logs/CollectMappedBasesPerSample.log"
    localrule: True
    shell:
        """
        cat {input} > {output}
        """

rule CatIdxStats_R:
    input:
        expand("idxstats/{sample}.idxstats.txt", sample=AllSamples)
    output:
        "QC/idxstats_combined.tsv"
    log:
        "logs/CatIdxStats_R.log"
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/CatIdxStats.R {output} {input} &> {log}
        """

rule QualimapRnaseq:
    input:
        gtf = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"),
        bam="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        bai="Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        results = "QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt",
        outdir = directory("QC/QualimapRnaseq/{sample}")
    log:
        "logs/QualimapRnaseq/{sample}.log"
    conda:
        "../envs/qualimap.yml"
    params:
        extra = "-p strand-specific-reverse"
    resources:
        mem_mb = 16000
    shell:
        """
        unset DISPLAY
        qualimap rnaseq -bam {input.bam} -gtf {input.gtf} {params.extra} --java-mem-size=12G -outdir {output.outdir}/ &> {log}
        """

rule MultiQC:
    input:
        expand("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam", sample= samples.index.unique()),
        expand("QC/QualimapRnaseq/{sample}/rnaseq_qc_results.txt", sample=samples.index.unique()),
        expand("FastqFastp/{sample}.fastp.json", sample=samples.index.unique())
    log: "logs/Multiqc.log"
    output:
        "Multiqc/multiqc_report.html"
    shell:
        """
        multiqc -f -o Multiqc/ Alignments/STAR_Align/ QC/QualimapRnaseq/ FastqFastp/ &> {log}
        """



rule CountReadsPerSample:
    input:
        expand("idxstats/{sample}.txt", sample=samples.index.unique())
    output:
        "../output/QC/ReadCountsPerSamples.tsv"
    log:
        "logs/CountReadsPerSample.log"
    shell:
        """
        # exec > {log} 2>&1
        # set -x
        for f in {input}
        do
           printf "%s\\t%s\\n" $f $(awk -F'\\t' '$1~"^chr[1-9]" {{sum+=$3}} END {{print sum}}' $f) >> {output}
        done
        """

rule CatIdxStats_R:
    input:
        expand("idxstats/{sample}.txt", sample=samples.index.unique()[0:10])
    output:
        "test_script.tsv"
    log:
        "logs/CatIdxStats_R.log"
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/CatIdxStats.R {output} {input} &> {log}
        """

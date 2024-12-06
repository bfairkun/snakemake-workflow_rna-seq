rule CopyFastq:
    """
    Useful for when a single sample is spread across multiple fastq, or when the original fastq is in long term storage cds
    """
    input:
        R1 = lambda wildcards: samples.loc[wildcards.sample]['R1'],
        R2 = lambda wildcards: samples.loc[wildcards.sample]['R2'],
    output:
        R1 = temp("Fastq/{sample}.R1.fastq.gz"),
        R2 = temp("Fastq/{sample}.R2.fastq.gz"),
    wildcard_constraints:
        sample = PairedEndSamples_wildcards_regex
    shell:
        """
        cat {input.R1} > {output.R1}
        cat {input.R2} > {output.R2}
        """

rule CopyFastq_SE:
    """
    Useful for when a single sample is spread across multiple fastq, or when the original fastq is in long term storage cds
    """
    input:
        R1 = lambda wildcards: samples.loc[wildcards.sample]['R1'],
    output:
        R1 = temp("Fastq/{sample}.R1.fastq.gz"),
    wildcard_constraints:
        sample = SingleEndSamples_wildcards_regex
    shell:
        """
        cat {input.R1} > {output.R1}
        """

rule fastp:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = "Fastq/{sample}.R1.fastq.gz",
        R2 = "Fastq/{sample}.R2.fastq.gz",
    output:
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        R2 = "FastqFastp/{sample}.R2.fastq.gz",
        html = "FastqFastp/{sample}.fastp.html",
        json = "FastqFastp/{sample}.fastp.json"
    params:
    resources:
        mem_mb = GetMemForSuccessiveAttempts(8000, 24000)
    wildcard_constraints:
        sample = PairedEndSamples_wildcards_regex
    log:
        "logs/fastp/{sample}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} -I {input.R2}  -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} &> {log}
        """

rule fastp_SE:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = "Fastq/{sample}.R1.fastq.gz",
    output:
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        html = "FastqFastp/{sample}.fastp.html",
        json = "FastqFastp/{sample}.fastp.json"
    params:
    resources:
        mem_mb = GetMemForSuccessiveAttempts(8000, 24000)
    wildcard_constraints:
        sample = SingleEndSamples_wildcards_regex
    log:
        "logs/fastp/{sample}.log"
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp -i {input.R1} -o {output.R1} --html {output.html} --json {output.json} &> {log}
        """

rule STAR_Align:
    input:
        index = lambda wildcards: config['GenomesPrefix'] + samples.loc[wildcards.sample]['STARGenomeName'] + "/STARIndex",
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        R2 = lambda wildcards: "FastqFastp/{sample}.R2.fastq.gz" if wildcards.sample in samples_PairedEnd else []
    output:
        outdir = directory("Alignments/STAR_Align/{sample}"),
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        align_log = "Alignments/STAR_Align/{sample}/Log.final.out"
    threads: 8
    log: "logs/STAR_Align/{sample}.log"
    params:
        GetSTARIndexDir = "/project2/yangili1/bjf79/ChromatinSplicingQTLs/code/ReferenceGenome/STARIndex/",
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
        extra = "--twopassMode Basic"
    resources:
        tasks = 9,
        mem_mb = 48000,
        # N = 1
    shell:
        """
        STAR --readMapNumber {params.readMapNumber} --outFileNamePrefix {output.outdir}/ --genomeDir {input.index}/ --readFilesIn {input.R1} {input.R2}  --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 --limitBAMsortRAM 8000000000 {params.ENCODE_params} --outSAMstrandField intronMotif {params.extra} &> {log}
        """


rule indexBam:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "logs/indexBam/{sample}.log"
    params:
        GetIndexingParamsFromSampleName
    output:
        index = touch("Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
    shell: "samtools index {params} {input} &> {log} && touch {output.index}"


rule idxstats:
    input:
        bam = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam",
        index = "Alignments/STAR_Align/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
    output:
        "idxstats/{sample}.idxstats.txt"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

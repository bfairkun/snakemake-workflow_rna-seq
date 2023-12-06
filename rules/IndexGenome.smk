rule DownloadFastaAndGtf:
    output:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
        # fa = "{{GenomesPrefix}}{GenomeName}/Reference.fa",
        # gtf = "{{GenomesPrefix}}{GenomeName}/Reference.gtf"
    params:
        fa_link = lambda wildcards: STAR_genomes.loc['GRCh38_GencodeRelease44Comprehensive']['FastaLink'],
        gtf_link = lambda wildcards: STAR_genomes.loc['GRCh38_GencodeRelease44Comprehensive']['GtfLink'],
    shell:
        """
        wget -O- {params.fa_link} | zcat > {output.fa}
        wget -O- {params.gtf_link} | zcat > {output.gtf}
        """

rule GetBasicGtf:
    """
    Some analysis are better done with just 'basic' tagged transcripts,
    transcript structures thought to encode functional transcripts. This shell
    command works with Gencode formatted gtf, haven't looked at Ensembl gtfs
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf"
    shell:
        """
        awk '/^##/ || /tag "basic"/ || $3=="gene"' {input.gtf} > {output.gtf}
        """


rule indexHg38Ref:
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
    output:
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    shell:
        """
        samtools faidx {input.fa}
        """

rule STAR_make_index:
    """
    did not work on bigmem2. Never figured out why (the log file didn't
    indicate anything). Ran on login node with success.
    """
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        index = directory(config['GenomesPrefix'] + "{GenomeName}/STARIndex")
    log:
        "logs/STAR_make_index/{GenomeName}.log"
    params:
    threads: 4
    resources:
        mem = "52000",
        # partition = "bigmem2",
        ntasks = 5
    shell:
        """
        mkdir -p {output.index}
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN {threads} --genomeDir {output.index}/ --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fa} &> {log}
        """

rule Gather_STAR_Indexes:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/STARIndex", GenomeName = STAR_genomes.index)

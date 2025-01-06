rule DownloadFastaAndGtf:
    output:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    params:
        fa_link = lambda wildcards: STAR_genomes.loc[wildcards.GenomeName]['FastaLink'],
        gtf_link = lambda wildcards: STAR_genomes.loc[wildcards.GenomeName]['GtfLink'],
    localrule: True
    shell:
        """
        wget -O- {params.fa_link} | zcat > {output.fa}
        wget -O- {params.gtf_link} | zcat > {output.gtf}
        """

rule GetBasicGtf:
    """
    Some analysis are better done with just 'basic' tagged transcripts,
    transcript structures thought to encode functional transcripts. I tried to
    make this shell command work for ensembl or Gencode gtfs, or maybe even
    RefSeq, but no guarantees. Also, remove any trailing whitespaces that might
    be in input gtf
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf"
    shell:
        """
        awk '/^#/ || /tag "basic"/ || /transcript_biotype "mRNA"/ || /transcript_biotype "protein_coding"/ || $3=="gene" || /transcript_id "[XN]M/' {input.gtf} | sed -e 's/[[:space:]]*$//' > {output.gtf}
        # check if output has any non comment lines, and just copy original gtf if not
        if ! grep -q -v -P '^#' {output.gtf}; then
            cat {input.gtf} | sed -e 's/[[:space:]]*$//' > {output.gtf}
        fi
        """

rule faidxGenome:
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
    output:
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    shell:
        """
        samtools faidx {input.fa}
        """

rule gtf_to_bed12:
    """
    This shell command works with Gencode formatted gtf, haven't looked at
    Ensembl gtfs
    """
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        bed = temp(config['GenomesPrefix'] + "{GenomeName}/Reference.bed.gz"),
        tbi = temp(config['GenomesPrefix'] + "{GenomeName}/Reference.bed.gz.tbi")
    conda:
        "../envs/bedparse.yml"
    shell:
        """
        bedparse gtf2bed {input.gtf} --extraFields gene_id,transcript_id,gene_type,gene_name,transcript_type,transcript_support_level,tag,transcript_name | awk -F'\\t' -v OFS='\\t' '{{$4=$NF; print $0}}' | bedtools sort -i - | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule colorize_gtfbed:
    """
    Only tested on human Gencode v44
    """
    input:
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.bed.gz",
    output:
        bed = temp(config['GenomesPrefix'] + "{GenomeName}/Reference.colored.bed"),
        colorkey = config['GenomesPrefix'] + "{GenomeName}/Reference.colored.ColorsKey.png"
    conda:
        "../envs/r_2.yml"
    log:
        "logs/colorize_gtfbed/{GenomeName}.log"
    shell:
        """
        Rscript scripts/ColorTranscriptsInGtfBed.R {input.bed} {output.bed} {output.colorkey} &> {log}
        """

rule index_colored_bed12:
    input:
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.colored.bed",
    output:
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.Transcripts.colored.bed.gz",
        tbi = config['GenomesPrefix'] + "{GenomeName}/Reference.Transcripts.colored.bed.gz.tbi"
    shell:
        """
        bedtools sort -i {input.bed} | bgzip /dev/stdin -c > {output.bed}
        tabix -p bed {output.bed}
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
        slurm_partition = "bigmem2",
        ntasks = 5
    shell:
        """
        mkdir -p {output.index}
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN {threads} --genomeDir {output.index}/ --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fa} &> {log}
        """

rule Gather_used_faidx:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai", GenomeName = samples['STARGenomeName'].unique())

rule Gather_STAR_Indexes:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/STARIndex", GenomeName = STAR_genomes.index)

rule Gather_used_STAR_Indexes:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/STARIndex", GenomeName = samples['STARGenomeName'].unique())

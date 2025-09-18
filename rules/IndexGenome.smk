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

rule SortIndexGtf:
    input:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf.gz",
        index  = touch(config['GenomesPrefix'] + "{GenomeName}/Reference.gtf.gz.indexing_done"),
    params:
        GetIndexingParamsFromGenomeName
    shell:
        """
        (grep '^#' {input.gtf} ; grep -v '^#' {input.gtf} | sort -k1,1 -k4,4n) | bgzip -c > {output.gtf} && tabix -f {params} -p gff {output.gtf} && touch {output.index}
        """

rule MakeIGV_GenomeFile:
    """
    Note that csi indexes don't seem to load with IGV genomes in my experience (indefenite loading). They can be opened as seperate tracks but not as part of the genome. Hence, the if statement in the jinja template
    """
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf.gz",
    output:
        config['GenomesPrefix'] + "{GenomeName}/Reference.igv.genome.json"
    params:
        name = lambda wildcards: STAR_genomes.loc[wildcards.GenomeName]['Notes'],
        index_suffix = GetIndexSuffix
    run:
        import jinja2
        # Read the template file
        template_content = """{
          "id": "{{name}}",
          "name": "{{name}}",
          "fastaURL": "Reference.fa",
          "indexURL": "Reference.fa.fai",
          "tracks": [
          {% if index_suffix == 'tbi' %}
            {
              "name": "Genes",
              "format": "gtf",
              "url": "Reference.gtf.gz",
              "indexURL": "Reference.gtf.gz.{{index_suffix}}",
              "indexed": true,
              "hidden" : false,
              "searchable": true,
              "removable": false
            }
          {% endif %}
          ]
        }
        """
        # Create a Jinja2 template from the content
        template = jinja2.Template(template_content)
        # Render the template with the wildcards and params
        rendered_content = template.render(id=wildcards.GenomeName, name=params.name, index_suffix = params.index_suffix)
        # Write the rendered content to the output file
        with open(output[0], 'w') as f:
            f.write(rendered_content)


rule MakeColored_Transcripts_bed:
    """
    Bed file for each transcript in gtf, colored by NMDFinderB status
    """
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.ColoredTranscripts.bed.gz",
        index  = touch(config['GenomesPrefix'] + "{GenomeName}/Reference.ColoredTranscripts.bed.gz.indexing_done"),
    params:
        tabixParams = GetIndexingParamsFromGenomeName
    shadow: "shallow"
    conda:
        "../scripts/leafcutter2/scripts/Reformat_gtf.conda_env.yml"
    shell:
        """
        python scripts/leafcutter2/scripts/Reformat_gtf.py -i {input.gtf} -fa {input.fa} -bed12 Transcripts.colored.bed  -o Transcripts.gtf -infer_gene_type_approach B -translation_approach C -gene_name_attribute_name gene_id
        bedtools sort -i Transcripts.colored.bed | bgzip -c /dev/stdin > {output.bed}
        tabix -f {params.tabixParams} -p bed {output.bed} && touch {output.index}
        """

rule Extract_introns:
    """
    Extract introns from gtf as bed file with 5'ss and 3'ss sequences
    """
    input:
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"
    output:
        bed = config['GenomesPrefix'] + "{GenomeName}/Reference.Introns.bed.gz",
        index  = touch(config['GenomesPrefix'] + "{GenomeName}/Reference.Introns.bed.gz.indexing_done"),
    log:
        "logs/Extract_introns/{GenomeName}.log"
    conda:
        "../envs/bio_python_tools.yml"
    params:
        tabixParams = GetIndexingParamsFromGenomeName
    shell:
        """
        (python scripts/ExtractIntronsFromGtf.py --gtf {input.gtf} --reference {input.fa} --output /dev/stdout |  awk -F'\\t' -v OFS='\\t' '{{split($4, a, "|"); print $1, $2, $3, a[2]"|"a[3], $5, $6}}' | sort | uniq | bedtools sort -i - | bgzip -c /dev/stdin > {output.bed} ) &> {log}
        tabix -f {params.tabixParams} -p bed {output.bed} && touch {output.index}
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
        expand(config['GenomesPrefix'] + "{GenomeName}/STARIndex", GenomeName = STAR_genomes.index),

rule Gather_ColoredTranscript_beds:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/Reference.ColoredTranscripts.bed.gz", GenomeName = STAR_genomes.index)

rule Gather_used_STAR_Indexes:
    input:
        expand(config['GenomesPrefix'] + "{GenomeName}/STARIndex", GenomeName = samples['STARGenomeName'].unique())

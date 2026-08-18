rule Make_new_samples_config:
    """
    Rewrite out samples.tsv sample config file, since pre-DAG calculation, links are determined from SRA accession if no local files or links to fastq files are given. Rather than looking up these links everytime the snakemake runs, we can just write a new sample config file with links already in place. Note that this also freezes the outcome of the bam prechecks: a bam that was found unusable will have been blanked in favor of its fastq/SRA fallback.
    """
    output:
        samples_tsv = "samples.SRA_accession_links_filled.tsv"
    run:
        samples.to_csv(output.samples_tsv, sep='\t', index=False)

rule StageProvidedBam:
    """
    Entry point for samples already aligned elsewhere. Verifies that chromosome names and lengths
    agree with the reference this sample will be quantified against, then symlinks the bam into
    place, or merges/sorts if it has to. rules/common.py performs the same check before the DAG is
    built when the reference is already indexed, but the first run against a new reference cannot,
    since Reference.fa.fai is itself built by this workflow.
    """
    input:
        bams = GetProvidedBams,
        fai = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai"),
    output:
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam"
    threads: 4
    resources:
        mem_mb = GetMemForSuccessiveAttempts(8000, 16000)
    log:
        "logs/StageProvidedBam/{sample}.log"
    wildcard_constraints:
        sample = wildcard_constraints_from_list(samples_from_bam)
    shell:
        """
        for bam in {input.bams}
        do
            samtools view -H $bam | awk -F'\\t' -v FAI={input.fai} -v BAM=$bam -v SAMPLE={wildcards.sample} '
                BEGIN {{ while ((getline line < FAI) > 0) {{ split(line, f, "\\t"); ref[f[1]] = f[2] }} }}
                $1=="@SQ" {{
                    name = ""; len = ""
                    for (i=2; i<=NF; i++) {{
                        if ($i ~ /^SN:/) name = substr($i, 4)
                        if ($i ~ /^LN:/) len = substr($i, 4)
                    }}
                    if (name in ref) {{ shared++; if (ref[name] != len) mismatched = mismatched " " name "(" len " vs " ref[name] ")" }}
                    else {{ absent = absent " " name }}
                }}
                END {{
                    if (shared == 0) {{
                        print "ERROR: sample " SAMPLE ": " BAM " has no chromosome names in common with " FAI ". Check for a chr-prefix naming mismatch. Fix the bam, or blank the bam column in the samples file to fall back to R1/R2/SRA_accession." > "/dev/stderr"
                        exit 1
                    }}
                    if (mismatched != "") {{
                        print "ERROR: sample " SAMPLE ": " BAM " chromosome lengths disagree with " FAI ":" mismatched ". This is likely a different assembly. Fix the bam, or blank the bam column in the samples file to fall back to R1/R2/SRA_accession." > "/dev/stderr"
                        exit 1
                    }}
                    if (absent != "") {{ print "Note: sample " SAMPLE ": contigs in bam absent from reference:" absent }}
                }}' >> {log}
        done

        if [ $(echo {input.bams} | wc -w) -eq 1 ]
        then
            if grep -q '^@HD.*SO:coordinate' <<< "$(samtools view -H {input.bams})"
            then
                ln -sr {input.bams} {output.bam} &>> {log}
            else
                samtools sort -@ {threads} -o {output.bam} {input.bams} &>> {log}
            fi
        else
            (samtools merge -u -@ {threads} - {input.bams} | samtools sort -@ {threads} -o {output.bam}) &>> {log}
        fi
        """

rule CopyFastq:
    """
    Useful for when a single sample is spread across multiple fastq, or when the original fastq is in long term storage cds
    """
    input:
        R1 = lambda wildcards: samples.loc[samples['sample']==wildcards.sample]['R1'],
        R2 = lambda wildcards: samples.loc[samples['sample']==wildcards.sample]['R2'],
    output:
        R1 = temp("Fastq/{sample}.R1.fastq.gz"),
        R2 = temp("Fastq/{sample}.R2.fastq.gz"),
    localrule: True
    wildcard_constraints:
        sample = wildcard_constraints_from_list(set(samples_from_local) & set(samples_PairedEnd))
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
        R1 = lambda wildcards: samples.loc[samples['sample']==wildcards.sample]['R1'],
    output:
        R1 = temp("Fastq/{sample}.R1.fastq.gz"),
    localrule: True
    wildcard_constraints:
        sample = wildcard_constraints_from_list(set(samples_from_local) & set(samples_SingleEnd))
    shell:
        """
        cat {input.R1} > {output.R1}
        """

rule DownloadFromAccession:
    """
    aspera_key is a param rather than an input so that it can be left blank in the config, in which
    case (or if the key file is missing) the fastq is fetched over ftp with wget instead
    """
    output:
        fastq = temp("Fastq/{sample}.{Read}.fastq.gz"),
    log:
        "logs/DownloadFastqFromAsperaLink/{sample}.{Read}.log"
    wildcard_constraints:
        Read = "R1|R2",
        sample = wildcard_constraints_from_list(samples_from_links)
    shadow: "shallow"
    conda:
        "../envs/ascp.yml"
    localrule: True
    params:
        aspera_key = config['aspera_key'],
        link = lambda wildcards: samples.loc[samples['sample']==wildcards.sample][wildcards.Read + '_link'].tolist(),
        aspera_link = lambda wildcards: [link.replace("ftp://ftp.sra.ebi.ac.uk/", "era-fasp@fasp.sra.ebi.ac.uk:") for link in samples.loc[samples['sample']==wildcards.sample][wildcards.Read + '_link'].tolist()]
    shell:
        """
        if [[ -e "{params.aspera_key}" && ! -z "{params.link}" ]]; then
            for link in {params.link}
            do
                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
                ascp -v -QT -l 300m -P33001 -i {params.aspera_key} {params.aspera_link} $tmpfile &>> {log}
                cat $tmpfile >> {output.fastq}
                rm $tmpfile
            done
        else
            for link in {params.link}
            do
                tmpfile=$(mktemp -p . tmp.download.XXXXXXXX.fastq.gz)
                wget -O $tmpfile ${{link}} &>> {log}
                cat $tmpfile >> {output.fastq}
                rm $tmpfile
            done
        fi
        """

rule fastp:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = "Fastq/{sample}.R1.fastq.gz",
        R2 = "Fastq/{sample}.R2.fastq.gz",
    output:
        R1 = temp("FastqFastp/{sample}.R1.fastq.gz"),
        R2 = temp("FastqFastp/{sample}.R2.fastq.gz"),
        html = "FastqFastp/{sample}.fastp.html",
        json = "FastqFastp/{sample}.fastp.json"
    params:
        cmd = lambda wildcards: "fastp" if samples.loc[samples['sample']==wildcards.sample]['Aligner'].tolist()[0]=='STAR' else "fastplong",
        extra = lambda wildcards: "" if samples.loc[samples['sample']==wildcards.sample]['Aligner'].tolist()[0]=='STAR' else "-Q"
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
        {params.cmd} {params.extra} -i {input.R1} -I {input.R2}  -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} &> {log}
        """

rule fastp_SE:
    """
    clips adapters, can handle UMIs
    """
    input:
        R1 = "Fastq/{sample}.R1.fastq.gz",
    output:
        R1 = temp("FastqFastp/{sample}.R1.fastq.gz"),
        html = "FastqFastp/{sample}.fastp.html",
        json = "FastqFastp/{sample}.fastp.json"
    params:
        cmd = lambda wildcards: "fastp" if samples.loc[samples['sample']==wildcards.sample]['Aligner'].tolist()[0]=='STAR' else "fastplong",
        extra = lambda wildcards: "" if samples.loc[samples['sample']==wildcards.sample]['Aligner'].tolist()[0]=='STAR' else "-Q"
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
        {params.cmd} {params.extra} -i {input.R1} -o {output.R1} --html {output.html} --json {output.json} &> {log}
        """

rule STAR_Align:
    input:
        index = lambda wildcards: config['GenomesPrefix'] + samples.loc[samples['sample']==wildcards.sample]['STARGenomeName'].tolist()[0] + "/STARIndex",
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        R2 = lambda wildcards: "FastqFastp/{sample}.R2.fastq.gz" if wildcards.sample in samples_PairedEnd.tolist() else []
    output:
        outdir = directory("Alignments/{sample}"),
        bam = temp("Alignments/{sample}/Aligned.out.bam"),
        align_log = "Alignments/{sample}/Log.final.out"
    threads: 8
    log: "logs/STAR_Align/{sample}.log"
    params:
        readMapNumber = -1,
        ENCODE_params = "--outFilterType BySJout --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
        extra = "--twopassMode Basic --outSAMattributes All --limitOutSJcollapsed 4000000"
    resources:
        mem_mb = 48000,
        # N = 1
    wildcard_constraints:
        sample = wildcard_constraints_from_list(samples_STAR_aligned_here)
    shell:
        """
        STAR --readMapNumber {params.readMapNumber} --outFileNamePrefix {output.outdir}/ --genomeDir {input.index}/ --readFilesIn {input.R1} {input.R2}  --outSAMtype BAM Unsorted --readFilesCommand zcat --runThreadN {threads} --outSAMmultNmax 1 {params.ENCODE_params} --outSAMstrandField intronMotif {params.extra} &> {log}
        """

rule sort_bam:
    input:
        bam = "Alignments/{sample}/Aligned.out.bam",
    output:
        sorted_bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam"
    threads: 6
    resources:
        mem_mb = 48000,
    wildcard_constraints:
        sample = wildcard_constraints_from_list(samples_needing_alignment)
    log:
        "logs/sort_bam/{sample}.log"
    shell:
        """
        samtools sort -@ {threads} -m 4000M -o {output.sorted_bam} {input.bam} &> {log}
        """

rule minimap2:
    input:
        R1 = "FastqFastp/{sample}.R1.fastq.gz",
        fa = lambda wildcards: config['GenomesPrefix'] + samples.loc[samples['sample']==wildcards.sample]['STARGenomeName'].tolist()[0] + "/Reference.fa",
    output:
        bam = temp("Alignments/{sample}/Aligned.out.bam"),
        align_log = "Alignments/{sample}/Log.final.out"
    resources:
        mem_mb = 48000,
    log: "logs/minimap2/{sample}.log"
    wildcard_constraints:
        sample = wildcard_constraints_from_list(samples.loc[(samples['Aligner']=='minimap2') & samples['bam'].isna()]['sample'])
    shadow: "shallow"
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 --secondary=no -a -x splice {input.fa} {input.R1} > Alignments.sam 2> {output.align_log}
        # Add strand tag (XS:A:[+-]) and remove antisense spliced reads (ts:A:-)
        (cat <(samtools view -H Alignments.sam) <(samtools view -F16 Alignments.sam | perl -lane 'if (!/\\tts:A:-/) {{print "$_\\tXS:A:+"}}') <(samtools view -f16 Alignments.sam | perl -lane 'if (!/\\tts:A:-/) {{print "$_\\tXS:A:-"}}') | samtools view -b -o {output.bam}) &>> {log}
        """

rule indexBam:
    input:
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "logs/indexBam/{sample}.log"
    params:
        GetIndexingParamsFromSampleName
    output:
        index = touch("Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done"),
    shell: "samtools index {params} {input} &> {log} && touch {output.index}"


rule idxstats:
    input:
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        index = "Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
    output:
        "idxstats/{sample}.idxstats.txt"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

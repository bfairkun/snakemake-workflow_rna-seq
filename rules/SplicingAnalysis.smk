rule ExtractJuncs:
    input:
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        index = "Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
    output:
        "SplicingAnalysis/juncfiles/{sample}.junc",
    params:
        strand = "0"
    conda:
        "../envs/regtools.yml"
    log:
        "logs/ExtractJuncs/{sample}.log"
    shell:
        """
        (regtools junctions extract -m 20 -s {params.strand} {input.bam} > {output}) &> {log}
        """


rule annotate_juncfiles:
    """
    this command breaks when gtf has trailing spaces, like in some NCBI gtfs
    https://github.com/griffithlab/regtools/issues/92
    """
    input:
        fa = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa"),
        fai = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai"),
        gtf = FillGenomeNameInFormattedString(config['GenomesPrefix'] + "{GenomeName}/Reference.gtf"),
        juncs = "SplicingAnalysis/juncfiles/{sample}.junc",
    output:
        counts = "SplicingAnalysis/juncfiles/{sample}.junccounts.tsv.gz"
    log:
        "logs/annotate_juncfiles/{sample}.log"
    conda:
        "../envs/regtools.yml"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 48000)
    shell:
        """
        (regtools junctions annotate {input.juncs} {input.fa} {input.gtf} | awk -F'\\t' -v OFS='\\t' 'NR>1 {{$4=$1"_"$2"_"$3"_"$6; print $4, $5}}' | gzip - > {output.counts} ) &> {log}
        """


rule ConcatJuncFilesAndKeepUniq:
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard("SplicingAnalysis/juncfiles/{sample}.junc"),
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc.gz"
    log:
        "logs/ConcatJuncFilesAndKeepUniq/{GenomeName}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 48000)
    shell:
        """
        (awk -v OFS='\\t' '{{ split($11, blockSizes, ","); JuncStart=$2+blockSizes[1]; JuncEnd=$3-blockSizes[2]; print $0, JuncStart, JuncEnd }}' {input} | sort -k1,1 -k6,6 -k13,13n -k14,14n | python scripts/AggregateSortedCattedJuncBeds.py | bedtools sort -i - | bgzip -c /dev/stdin  > {output}) &> {log}
        """

rule AnnotateConcatedUniqJuncFile_basic:
    """
    note that regtools end coordinate is off by one. The right coordinate should be adjusted down 1 for proper viewing of junc in IGV
    """
    input:
        junc = "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.junc.gz",
        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.basic.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa"
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz"
    log:
        "logs/AnnotateConcatedUniqJuncFile_hg38Basic.{GenomeName}.log"
    conda:
        "../envs/regtools.yml"
    shell:
        """
        (regtools junctions annotate {input.junc} {input.fa} {input.gtf} | gzip - > {output} ) &> {log}
        """

rule Add_splice_site_scores_to_regtools_annotate:
    input:
        annotated_junctions="SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz",
        reference_fasta=config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
        AnnotatedIntronsWithSS = config['GenomesPrefix'] + "{GenomeName}/Reference.Introns.bed.gz"
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.with_ss_scores.tsv.gz"
    conda:
        "../scripts/leafcutter2/scripts/Reformat_gtf.conda_env.yml"
    log:
        "logs/Add_splice_site_scores_to_regtools_annotate.{GenomeName}.log"
    shell:
        """
        python scripts/Add_SS_To_RegtoolsAnnotate.py \
            --input {input.annotated_junctions} \
            --reference {input.reference_fasta} \
            --introns {input.AnnotatedIntronsWithSS} \
            --output {output} &> {log}
        """

rule make_leafcutter_juncfile:
    input:
        ExpandAllSamplesInFormatStringFromGenomeNameWildcard("SplicingAnalysis/juncfiles/{sample}.junc"),
    output:
        "SplicingAnalysis/leafcutter/{GenomeName}/juncfilelist.txt"
    params:
        SamplesToRemove = ""
    run:
        import os
        if params.SamplesToRemove:
            SamplesToRemove = open(params.SamplesToRemove, 'r').read().split('\n')
        else:
            SamplesToRemove=[]
        with open(output[0], "w") as out:
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename not in  SamplesToRemove:
                    out.write(filepath + '\n')

rule leafcutter_cluster:
    input:
        juncs = ExpandAllSamplesInFormatStringFromGenomeNameWildcard("SplicingAnalysis/juncfiles/{sample}.junc"),
        juncfile_list = "SplicingAnalysis/leafcutter/{GenomeName}/juncfilelist.txt"
    output:
        outdir = directory("SplicingAnalysis/leafcutter/{GenomeName}/clustering/"),
        counts = "SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind.counts.gz",
        numers = "SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind_numers.counts.gz"
    shadow: "shallow"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 48000)
    log:
        "logs/leafcutter_cluster/{GenomeName}.log"
    params:
        "-p 0.0001"
    shell:
        """
        python scripts/leafcutter/clustering/leafcutter_cluster_regtools.py -j {input.juncfile_list} {params} -r {output.outdir} -k True &> {log}
        """

rule leafcutter_to_PSI:
    input:
        numers = "SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind_numers.counts.gz"
    output:
        juncs = temp("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/JuncCounts.bed"),
        PSI = temp("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/PSI.bed"),
        PSIByMax = temp("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/PSI_ByMax.bed"),
        PSIDenom = temp("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/PSIDenom.bed"),
        PSIByMaxDenom = temp("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/PSI_ByMaxDenom.bed")
    log:
        "logs/leafcutter_to_PSI/{GenomeName}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 54000)
    conda:
        "../envs/r_2.yml"
    shell:
        """
        Rscript scripts/leafcutter_to_PSI.R {input.numers} {output.PSI} {output.PSIByMax} {output.juncs} {output.PSIDenom} {output.PSIByMaxDenom} &> {log}
        """

rule bgzip_PSI_bed:
    input:
        bed = "SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/{Metric}.bed",
    output:
        bed = "SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/{Metric}.sorted.bed.gz",
        index  = touch("SplicingAnalysis/leafcutter/{GenomeName}/juncTableBeds/{Metric}.sorted.bed.gz.indexing_done"),
    log:
        "logs/bgzip_PSI_bed/{GenomeName}/{Metric}.log"
    resources:
        mem_mb = GetMemForSuccessiveAttempts(24000, 54000)
    params:
        GetIndexingParamsFromGenomeName
    shell:
        """
        (bedtools sort -header -i {input.bed} | bgzip /dev/stdin -c > {output.bed}) &> {log}
        (tabix {params} -f -p bed {output.bed}) &>> {log} && touch {output.index}
        """

rule Get5ssSeqs:
    """
    Filtered out entries with N in sequence
    """
    input:
        basic = "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.tsv.gz",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
    output:
        "SplicingAnalysis/ObservedJuncsAnnotations/{GenomeName}.uniq.annotated.DonorSeq.tsv"
    shell:
        """
        zcat {input.basic} | awk -v OFS='\\t' -F'\\t' 'NR>1 {{print $1, $2, $3, $1"_"$2"_"$3"_"$6, ".", $6}}' | sort -u | awk -v OFS='\\t' -F'\\t'  '$6=="+" {{$2=$2-4; $3=$2+11; print $0}} $6=="-" {{$3=$3+3; $2=$3-11; print $0}}' | bedtools getfasta -tab -bed - -s -name -fi {input.fa} | grep -v 'N' > {output}
        """

rule SpliceQ:
    input:
        bam = "Alignments/{sample}/Aligned.sortedByCoord.out.bam",
        index = "Alignments/{sample}/Aligned.sortedByCoord.out.bam.indexing_done",
        gtf = lambda wildcards: config['GenomesPrefix'] + samples.loc[wildcards.sample]['STARGenomeName'] + "/Reference.basic.gtf",
    output:
        "SplicingAnalysis/SpliceQ/{sample}.tsv.gz"
    log:
        "logs/SpliceQ/{sample}.log"
    params:
        "-c 0 -i"
    threads:
        4
    shadow:
        "shallow"
    conda:
        "../envs/spliceq.yml"
    shell:
        """
        exec 2> {log}
        tmpfile=$(mktemp -p . tmp.{wildcards.sample}.XXXXXXXX.txt)
        SPLICE-q.py {params} -b {input.bam} -g {input.gtf} -p {threads} -o $tmpfile
        gzip -c $tmpfile > {output}
        """

rule ConvertLeacfutterJuncNames_To_True_Bed_Coords:
    input:
        junclist = "SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind.counts.gz",
    output:
        "SplicingAnalysis/ClassifyJuncs/{GenomeName}.juncList.tsv.gz"
    shell:
        """
        zcat {input} | awk 'NR==1 {{ print $1 }} NR>1 {{split($1, a, ":"); print a[1]":"a[2]":"a[3]-1":"a[4] }}' | gzip - > {output}
        """


rule leafcutter2_ClassifyJuncs_ClusterPerInd:
    """
    Works with Gencode GTFs only. Need to reformat GTF for others
    """
    input:
        junclist = "SplicingAnalysis/ClassifyJuncs/{GenomeName}.juncList.tsv.gz",

        gtf = config['GenomesPrefix'] + "{GenomeName}/Reference.gtf",
        fa = config['GenomesPrefix'] + "{GenomeName}/Reference.fa",
        fai = config['GenomesPrefix'] + "{GenomeName}/Reference.fa.fai",
    output:
        outdir = directory("SplicingAnalysis/ClassifyJuncs/{GenomeName}"),
        Classifications = "SplicingAnalysis/ClassifyJuncs/{GenomeName}/Leaf2_junction_classifications.txt",
    log:
        "logs/leafcutter2_ClassifyJuncs_ClusterPerInd/{GenomeName}.log"
    params:
        rundir = "MazinLeafcutterAnalysis/ClassifyJuncs"
    conda:
        "../envs/bedparse.yml"
    resources:
        mem_mb = 16000
    shell:
        """
        python scripts/leafcutter2/scripts/SpliceJunctionClassifier.py -c {input.junclist} -G {input.fa} -A {input.gtf} -v -r {output.outdir} &> {log}
        """

rule leafcutter_ds_contrasts:
    input:
        groupfile = lambda wildcards: os.path.abspath(config['contrast_group_files_prefix'] + "{contrast}.txt"),
        numers = lambda wildcards: get_filled_path_from_contrast(wildcards, "SplicingAnalysis/leafcutter/{GenomeName}/clustering/leafcutter_perind_numers.counts.gz"),
    output:
        directory("SplicingAnalysis/differential_splicing/{contrast}/")
    threads: 4
    wildcard_constraints:
        treatment = "|".join(contrasts)
    resources:
        ntasks = 5,
        mem_mb = 24000
    params:
        Prefix = "MazinLeafcutterAnalysis/Contrasts_ds/",
        ExtraParams = "-i 2 -g 2"
    log:
        "logs/leafcutter_ds/{contrast}.log"
    shell:
        """
        mkdir -p {output}
        /software/R-3.4.3-el7-x86_64/bin/Rscript scripts/leafcutter/scripts/leafcutter_ds.R -p {threads} -o {output}/leaf {params.ExtraParams} {input.numers} {input.groupfile} &> {log}
        """
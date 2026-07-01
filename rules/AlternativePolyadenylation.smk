# Alternative polyadenylation (APA) / 3'UTR isoform quantification via QAPA + Salmon.
#
# Config keys read from config['apa_qapa'] (nest this under your rna_seq config):
#
#   apa_qapa:
#     genome: GRCh38_GencodeRelease44Comprehensive   # genome name under GenomesPrefix
#     scratch: "/scratch/midway3/bjf79/"              # large Salmon outputs go here
#     extra_polya_bed: ""                             # optional custom polyA BED
#     filter_gtf_tags: ["mRNA_end_NF"]                # transcript tags to exclude; [] = no filter
#     groups:                                         # one PAU matrix per group
#       MyGroup:
#         Treatment: ["Branaplam", "DMSO"]
#
# Sample groups are defined as column→values filters on your samples.tsv, OR
# as an explicit 'samples' list.  The final output per group is:
#   QAPA/{group}/PAU_matrix.tsv.gz

wildcard_constraints:
    apa_group = "[^/]+"


def apa_samples(group):
    """Return sample IDs for a given apa_qapa group."""
    cfg = config["apa_qapa"]["groups"][group]
    if "samples" in cfg:
        return list(cfg["samples"])
    result = samples.copy()
    for col, vals in cfg.items():
        result = result[result[col].astype(str).isin([str(v) for v in vals])]
    return result.index.unique().tolist()


_APA_GENOME  = config["apa_qapa"]["genome"]
_APA_REF_DIR = config["GenomesPrefix"] + _APA_GENOME + "/"
_APA_SCRATCH = config["apa_qapa"]["scratch"]
_APA_ANNO    = "QAPA/Annotation/" + _APA_GENOME


# ── Reference build (once per genome) ─────────────────────────────────────────

rule QAPA_FilterGTF:
    """Drop GTF lines matching any tag in config['apa_qapa']['filter_gtf_tags'].
    Set to [] to pass through unchanged. Multiple tags joined as egrep alternation.
    """
    input:
        gtf = _APA_REF_DIR + "Reference.gtf"
    output:
        _APA_ANNO + ".filtered.gtf"
    params:
        pattern = lambda wc: "|".join(config["apa_qapa"].get("filter_gtf_tags", []))
    log:
        "logs/QAPA/FilterGTF.log"
    shell:
        """
        if [ -z "{params.pattern}" ]; then
            ln -sf {input.gtf} {output} 2> {log}
        else
            grep -vE "{params.pattern}" {input.gtf} > {output} 2> {log}
        fi
        """


rule QAPA_Identifiers:
    """Parse filtered GTF → QAPA gene/transcript identifiers table."""
    input:
        gtf = _APA_ANNO + ".filtered.gtf"
    output:
        _APA_ANNO + ".identifiers.txt"
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/Identifiers.log"
    shell:
        "python scripts/qapa_build_identifiers.py {input.gtf} {output} &> {log}"


rule QAPA_GenePred:
    """Convert filtered GTF → genePred format required by qapa build."""
    input:
        gtf = _APA_ANNO + ".filtered.gtf"
    output:
        _APA_ANNO + ".genePred"
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/GenePred.log"
    shell:
        "gtfToGenePred -genePredExt {input.gtf} {output} 2> {log}"


rule QAPA_Build:
    """Build 3'UTR isoform BED from annotation (GTF-only, or with extra polyA sites)."""
    input:
        genepred = _APA_ANNO + ".genePred",
        db       = _APA_ANNO + ".identifiers.txt"
    output:
        _APA_ANNO + ".3UTRs.bed"
    params:
        polya_flag = lambda wc: (
            "-o " + config["apa_qapa"]["extra_polya_bed"]
            if config["apa_qapa"].get("extra_polya_bed", "")
            else "-N"
        )
    resources:
        mem_mb = 16000
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/Build.log"
    shell:
        "qapa build {params.polya_flag} --db {input.db} {input.genepred} > {output} 2> {log}"


rule QAPA_Fasta:
    """Extract 3'UTR isoform sequences for Salmon indexing."""
    input:
        bed = _APA_ANNO + ".3UTRs.bed",
        fa  = _APA_REF_DIR + "Reference.fa"
    output:
        _APA_ANNO + ".3UTRs.fa"
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/Fasta.log"
    shell:
        "qapa fasta -f {input.fa} {input.bed} {output} &> {log}"


rule QAPA_SalmonIndex:
    """Build Salmon index over 3'UTR isoform sequences."""
    input:
        _APA_ANNO + ".3UTRs.fa"
    output:
        directory(_APA_ANNO + ".salmon_index")
    threads: 8
    resources:
        mem_mb = 8000
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/SalmonIndex.log"
    shell:
        "salmon index -t {input} -i {output} -p {threads} &> {log}"


# ── Per-sample Salmon quantification ──────────────────────────────────────────

rule QAPA_SalmonQuant:
    """Align-free Salmon quant of 3'UTR isoforms from existing FASTQs.

    Uses ancient() so upstream alignment rules are never re-triggered.
    Outputs land in scratch (not git-tracked); clean up manually when done.
    Library type auto-detected (-l A).
    """
    input:
        r1    = ancient("Fastq/{sample}.R1.fastq.gz"),
        r2    = ancient("Fastq/{sample}.R2.fastq.gz"),
        index = _APA_ANNO + ".salmon_index"
    output:
        _APA_SCRATCH + "QAPA/salmon/{sample}/quant.sf"
    threads: 8
    resources:
        mem_mb = 8000
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/SalmonQuant/{sample}.log"
    shell:
        """
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} \
            -o $(dirname {output}) \
            &> {log}
        """


# ── Per-group PAU aggregation ──────────────────────────────────────────────────

rule QAPA_Quant:
    """Aggregate per-sample Salmon results into a PAU matrix for a sample group.

    Output columns: APA_ID, Transcript, Gene, Gene_Name, genomic coords,
    then <sample>.PAU and <sample>.TPM per sample.
    PAU (Poly-A Usage) = relative isoform proportion within each gene [0–100].
    Memory scales with samples: compute_pau.R dcasts large matrices for
    large groups; 16 GB covers 500+ samples.
    """
    input:
        quant_sfs = lambda wc: expand(
            _APA_SCRATCH + "QAPA/salmon/{sample}/quant.sf",
            sample=apa_samples(wc.apa_group)
        ),
        db = _APA_ANNO + ".identifiers.txt"
    output:
        "QAPA/{apa_group}/PAU_matrix.tsv.gz"
    resources:
        mem_mb = 16000
    conda:
        "../envs/qapa.yaml"
    log:
        "logs/QAPA/Quant/{apa_group}.log"
    shell:
        "qapa quant --db {input.db} {input.quant_sfs} 2> {log} | gzip > {output}"


rule QAPA_all:
    """Collect PAU matrices for all configured apa_qapa groups."""
    input:
        expand(
            "QAPA/{apa_group}/PAU_matrix.tsv.gz",
            apa_group=list(config["apa_qapa"]["groups"].keys())
        )

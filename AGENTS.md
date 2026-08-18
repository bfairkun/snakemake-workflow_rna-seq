# Repository conventions

Snakemake workflow for RNA-seq alignment and quantification, usually consumed as a submodule of a project workflow. See `README.md` for configuration and usage.

## Layout

* `Snakefile` — config, wildcard constraints, includes, and `rule all`.
* `rules/common.py` — samples/genomes tables, entry-point resolution, sample subsets, and input functions. Included first; everything else depends on the names it defines.
* `rules/*.smk` — rules grouped by stage. `SplicingAnalysis.smk` inherits rules from `ExpressionAnalysis.smk`, so include order matters.
* `envs/*.yml` — one conda env per tool, referenced by `conda:` in the rule that needs it.
* `scripts/` — helper scripts plus submodules (leafcutter, leafcutter2, GenometracksByGenotype). Run `git submodule update --init --recursive` after cloning.

## Conventions

* Keep `rules/common.py` a flat script of pandas operations and plain functions. Sample subsets are computed once at parse time and used as wildcard constraints via `wildcard_constraints_from_list`.
* Rules that fan a sample list into one job take an input function from the `ExpandAll*FromGenomeName*` family rather than expanding inline, since samples are grouped by genome and strandedness.
* Every rule writes to `logs/{RuleName}/{wildcards}.log`. Errors a user must act on go to stderr, not just the log, so they surface in the Snakemake output.
* Memory for rules that sometimes need more on a retry comes from `GetMemForSuccessiveAttempts(first, second, ...)`.
* Comment sparingly. Rules carry a short docstring only when the name doesn't convey the purpose; helper functions in `common.py` generally carry none.

## Samples file contract

`config/samples.tsv` and `config/STAR_Genome_List.tsv` are validated against `schemas/`. Adding a column means adding it to the schema with a description — the schema is the reference documentation for those files, and `snakemake.utils.validate` also fills column defaults from it.

Sample names must satisfy `^[A-Za-z_][A-Za-z0-9_.]*$`. Hyphens are excluded deliberately: sample names become count-table column names read by R, which mangles them.

Rows sharing a `sample` value are combined (fastq concatenated, bams merged). Entry-point precedence is `bam` > `R1`/`R2` > `SRA_accession`; a sample must not mix a bam row with a fastq row.

## Testing changes

    ./.test/make_test_data.sh
    snakemake -n --configfile .test/config.yaml

Covers all three entry points. A change touching entry-point resolution should also be checked against the failure paths: unusable bam with and without a fastq fallback, and a bam-entry row missing `Library_Layout`.

## Cluster constraints (RCC Midway3)

Measured, not assumed — these are why some rules are pinned and others are not:

* Compute nodes reach **port 443 but not port 21**. ENA reports `ftp://` links, so `GetDownloadLinks` rewrites them to `https://` (same host and path) and `DownloadFromAccession` runs on the cluster. Leaving them as ftp makes the job hang until it times out.
* **`/cds` is not mounted on compute nodes.** `CopyFastq` and `CopyFastq_SE` read from there, so they must stay `localrule: True`.
* A single connection to EBI sustains **~1.2 MiB/s** regardless of protocol (https 1.19, ftp 1.20) or node type. Aggregate scales linearly with concurrency and had not saturated at 24 streams (33.7 MiB/s). Downloads are therefore parallelism-bound, not bandwidth-bound.
* `DownloadFastaAndGtf` is left `localrule: True` deliberately: it is a handful of jobs per genome, and one link in `STAR_Genome_List.tsv` is plain `http://`, whose port has not been checked from a compute node.
* Aspera is not used; `era-fasp@fasp.sra.ebi.ac.uk` fails to authenticate with the standard public key (ascp 3.9.1), and parallel https beats it anyway.

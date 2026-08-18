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

## Known issues

* `rule DownloadFromAccession` is wrong for a sample whose reads are spread over several runs *when downloading via aspera*. The shell loops `for link in {params.link}` but then passes `{params.aspera_link}` — the whole list — to a single `ascp` call, so every iteration re-fetches all links into the same temp file. The wget branch below it handles the same case correctly, one link per iteration. The fix is to transform the loop variable in bash (`${link/ftp:\/\/ftp.sra.ebi.ac.uk\//era-fasp@fasp.sra.ebi.ac.uk:}`) and drop the now-unused `aspera_link` param; it needs an aspera key plus a multi-run accession to test against.

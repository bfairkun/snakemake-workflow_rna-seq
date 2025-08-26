rule generate_GenomeTracks_group_files:
    output:
        "Plot_contrasts/all_contrasts.groups.tsv"
    run:
        outdir = os.path.dirname(output[0])
        os.makedirs(outdir, exist_ok=True)
        output_prefix = output[0].split("Plot_contrasts")[0]
        rows = []
        for contrast in contrasts:
            df = contrasts_df[contrasts_df["ContrastName"] == contrast].copy()
            group_order = df["Group"].unique()
            color_map = {group_order[0]: "#bdbdbd"}
            if len(group_order) > 1:
                color_map[group_order[1]] = "#e41a1c"
            for group in group_order:
                group_df = df[df["Group"] == group]
                # Use the first STARGenomeName for the group (assuming all are the same)
                star_genome = group_df["STARGenomeName"].iloc[0]
                group_label = f"{group}.{contrast}"
                color = color_map[group]
                bedgz = os.path.join(output_prefix, f"SplicingAnalysis/leafcutter/{star_genome}/juncTableBeds/PSI.sorted.bed.gz")
                supergroup = contrast
                rows.append([group_label, color, bedgz, supergroup])
        out_df = pd.DataFrame(rows, columns=["Group_label", "color", "Bedgz", "Supergroup"])
        out_df.to_csv(output[0], sep="\t", index=False)

rule generate_all_samples_bigwig_table:
    output:
        "Plot_contrasts/all_contrasts_all_samples_bigwig_table.tsv"
    run:
        # Get the output prefix (everything up to 'Plot_contrasts')
        output_prefix = output[0].split("Plot_contrasts")[0]

        df = contrasts_df.copy()
        df["sample"] = df.apply(lambda row: f"{row['sample']}.{row['Group']}.{row['ContrastName']}", axis=1)
        df["bigwigPath"] = df.apply(
            lambda row: os.path.join(output_prefix, f"bigwigs/unstranded/{row['sample'].split('.')[0]}.bw"), axis=1
        )
        df["Group_label"] = df.apply(lambda row: f"{row['Group']}.{row['ContrastName']}", axis=1)
        df["strand"] = "."
        out_df = df[["sample", "bigwigPath", "Group_label", "strand"]]
        outdir = os.path.dirname(output[0])
        os.makedirs(outdir, exist_ok=True)
        out_df.to_csv(output[0], sep="\t", index=False)
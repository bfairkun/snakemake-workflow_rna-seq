rule generate_groups_files:
    output:
        expand("config/differential_splicing_expression_group_files/{contrast}.txt", contrast=contrast_names)
    run:
        for contrast in contrast_names:
            contrast_df = contrasts_df[contrasts_df["ContrastName"] == contrast]
            contrast_df[["sample", "Group"]].to_csv(f"config/differential_splicing_expression_group_files/{contrast}.txt", sep="\t", index=False, header=False)

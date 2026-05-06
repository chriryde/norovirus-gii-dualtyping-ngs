    # PARTIAL_RDRP_FILTERED_FASTA: Path = ct.extract_cds_feature_to_fasta(
    #     PARTIAL_FASTA, 
    #     partial_rdrp_filtered_sequences,
    #     PROJECT_ROOT / Path(config["paths"]["filtered_rdrp_fasta"]),
    #     "rdrp"
    # )

    # PARTIAL_RDRP_FILTERED_FASTA = PROJECT_ROOT / Path(PARTIAL_RDRP_FILTERED_FASTA)

    # PARTIAL_VP1_FILTERED_FASTA: Path = ct.extract_cds_feature_to_fasta(
    #     PARTIAL_FASTA, 
    #     partial_vp1_filtered_sequences,
    #     PROJECT_ROOT / Path(config["paths"]["filtered_vp1_fasta"]),
    #     "vp1"
    # )

    # PARTIAL_VP1_FILTERED_FASTA = PROJECT_ROOT / Path(PARTIAL_VP1_FILTERED_FASTA)

    # PARTIAL_JUNCTION_FILTERED_FASTA: Path = ct.extract_junction_sequences_to_fasta(
    #     PARTIAL_FASTA,
    #     PARTIAL_REGION_INFORMATION,
    #     partial_junction_filtered_sequences,
    #     PROJECT_ROOT / Path(config["paths"]["filtered_junction_fasta"])
    # )

    # PARTIAL_JUNCTION_FILTERED_FASTA = PROJECT_ROOT / Path(PARTIAL_JUNCTION_FILTERED_FASTA)

    # ct.typing_tool_intialise(
    #     PARTIAL_RDRP_FILTERED_FASTA,
    #     "partial_rdrp"
    # )

    # ct.typing_tool_intialise(
    #     PARTIAL_VP1_FILTERED_FASTA,
    #     "partial_vp1"
    # )

    # ct.typing_tool_intialise(
    #     PARTIAL_JUNCTION_FILTERED_FASTA,
    #     "partial_junction"
    # )

# PARTIAL_RDRP_FILTERED_FASTA = PROJECT_ROOT / Path(config["paths"]["filtered_rdrp_fasta"])
# PARTIAL_RDRP_TYPING_CSV: Path | None = ct.typing_tool_get_results(
#     TEMP_TYPING_DIR,
#     "partial_rdrp"
# )
# if PARTIAL_RDRP_TYPING_CSV is None:
#     print(f"Typing tool results are not finished yet. Halting pipeline.")
#     exit

# PARTIAL_RDRP_TYPING_CSV = PROJECT_ROOT / Path(PARTIAL_RDRP_TYPING_CSV)

# PARTIAL_VP1_FILTERED_FASTA = PROJECT_ROOT / Path(config["paths"]["filtered_vp1_fasta"])
# PARTIAL_VP1_TYPING_CSV: Path | None = ct.typing_tool_get_results(
#     TEMP_TYPING_DIR,
#     "partial_vp1"
# )
# if PARTIAL_VP1_TYPING_CSV is None:
#     print(f"Typing tool results are not finished yet. Halting pipeline.")
#     exit

# PARTIAL_VP1_TYPING_CSV = PROJECT_ROOT / Path(PARTIAL_VP1_TYPING_CSV)

# PARTIAL_JUNCTION_FILTERED_FASTA = PROJECT_ROOT / Path(config["paths"]["filtered_junction_fasta"])
# PARTIAL_JUNCTION_TYPING_CSV: Path | None = ct.typing_tool_get_results(
#     TEMP_TYPING_DIR,
#     "partial_junction"
# )
# if PARTIAL_JUNCTION_TYPING_CSV is None:
#     print(f"Typing tool results are not finished yet. Halting pipeline.")
#     exit

# PARTIAL_JUNCTION_TYPING_CSV = PROJECT_ROOT / Path(PARTIAL_JUNCTION_TYPING_CSV)

# partial_rdrp_filtered_fasta: data/01_curated/fasta/partial_rdrp_filtered.fna
# partial_vp1_filtered_fasta: data/01_curated/fasta/partial_vp1_filtered.fna
# partial_junction_filtered_fasta: data/01_curated/fasta/partial_junction_filtered.fna


# partial_rdrp_full_specification: data/01_curated/metadata/partial_rdrp_full_specification.csv
# partial_vp1_full_specification: data/01_curated/metadata/partial_vp1_full_specification.csv
# partial_junction_full_specification: data/01_curated/metadata/partial_junction_full_specification.csv
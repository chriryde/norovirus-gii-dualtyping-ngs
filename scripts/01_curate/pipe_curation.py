import curation_tools as ct
import yaml

from pathlib import Path

def project_path(pathlike: str | Path) -> Path:
    path = Path(pathlike)
    return path if path.is_absolute() else PROJECT_ROOT / path

CURATE_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (CURATE_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

INPUT_METADATA = project_path(config["paths"]["input_metadata"])
RAW_FASTA = project_path(config["paths"]["input_fasta"])
INPUT_CDS_FASTA = project_path(config["paths"]["input_cds_fasta"])
TEMP_TYPING_DIR = project_path(config["paths"]["temp_typing_dir"])

COMPLETE_FASTA_PATH = project_path(config["paths"]["complete_fasta"])
COMPLETE_FILTERED_FASTA_PATH = project_path(config["paths"]["complete_filtered_fasta"])
COMPLETE_CDS_FASTA_PATH = project_path(config["paths"]["complete_cds_fasta"])

FINAL_COMPLETE_FASTA = project_path(config["paths"]["final_complete_fasta"])

PARTIAL_FASTA_PATH = project_path(config["paths"]["partial_fasta"])
PARTIAL_FILTERED_FASTA_PATH = project_path(config["paths"]["partial_filtered_fasta"])
PARTIAL_CDS_FASTA_PATH = project_path(config["paths"]["partial_cds_fasta"])

FINAL_PARTIAL_FASTA = project_path(config["paths"]["final_partial_fasta"])

if TEMP_TYPING_DIR.is_dir():
    print("shouldn't run this")

    complete_sequences, partial_sequences = ct.get_sequence_completeness(
        INPUT_METADATA
    )

#|===| START OF COMPLETE BLOCK |===============================================|
    COMPLETE_FASTA: Path = ct.extract_and_save_to_fasta(
        RAW_FASTA, 
        complete_sequences,
        COMPLETE_FASTA_PATH
    )

    complete_gii_sequences: set[str] = ct.get_GII_sequences(
        INPUT_METADATA,
        complete_sequences
    )
    assert(complete_gii_sequences)

    complete_length_filtered_sequences: set[str] = ct.get_length_filtered_sequences(
        COMPLETE_FASTA,
        config["criteria"]["length_threshold"]
    )
    assert(complete_length_filtered_sequences)

    complete_unique_sequences: set[str] = ct.get_unique_sequences(
        COMPLETE_FASTA
    )
    assert(complete_unique_sequences)

    complete_annotated_sequences: set[str] = ct.get_annotated_sequences(
        INPUT_METADATA
    )
    assert(complete_annotated_sequences)

    complete_non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        COMPLETE_FASTA,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(complete_non_ambiguous_sequences)

    complete_filtered_sequences: set[str] = complete_sequences \
        & complete_length_filtered_sequences & complete_unique_sequences \
        & complete_annotated_sequences & complete_non_ambiguous_sequences \
        & complete_gii_sequences
    assert complete_filtered_sequences

    COMPLETE_FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        COMPLETE_FASTA, 
        complete_filtered_sequences,
        COMPLETE_FILTERED_FASTA_PATH
    )

    excluded_accessions = complete_sequences - complete_filtered_sequences
    ct.export_excluded_sequences(
        project_path(config["paths"]["excluded_csv"]),
        excluded_accessions=excluded_accessions,
        gii_sequences=complete_gii_sequences,
        complete_sequences=complete_sequences,
        length_filtered_sequences=complete_length_filtered_sequences,
        unique_sequences=complete_unique_sequences,
        annotated_sequences=complete_annotated_sequences,
        non_ambiguous_sequences=complete_non_ambiguous_sequences
    )

    # ct.typing_tool_intialise(
    #     COMPLETE_FILTERED_FASTA,
    #     "complete"
    # )

    COMPLETE_CDS_FASTA: Path = ct.extract_and_save_to_fasta(
        INPUT_CDS_FASTA, 
        complete_filtered_sequences,
        COMPLETE_CDS_FASTA_PATH
    )
#|===| END OF COMPLETE BLOCK |================================================|
#|                                                                            |
#|===| START OF PARTIAL BLOCK |===============================================|
    PARTIAL_FASTA: Path = ct.extract_and_save_to_fasta(
        RAW_FASTA, 
        partial_sequences,
        PARTIAL_FASTA_PATH
    )

    PARTIAL_CDS_FASTA: Path = ct.extract_and_save_to_fasta(
        INPUT_CDS_FASTA, 
        partial_sequences,
        PARTIAL_CDS_FASTA_PATH
    )

    PARTIAL_REGION_INFORMATION: Path | None = ct.get_genomic_info(
        project_path(config["paths"]["partial_region_information"]),
        typing_information=False,
        region_information=True,
        CDS_FASTA=PARTIAL_CDS_FASTA
    )
    if PARTIAL_REGION_INFORMATION is None:
        raise RuntimeError("Failed to create partial region information")

    partials_rdrp: set[str] = ct.get_rdrp_sequences(
        PARTIAL_REGION_INFORMATION,
    )

    partials_vp1: set[str] = ct.get_vp1_sequences(
        PARTIAL_REGION_INFORMATION,
    )

    partials_junction: set[str] = ct.get_junction_sequences(
        PARTIAL_REGION_INFORMATION,
    )

    partial_gii_sequences: set[str] = ct.get_GII_sequences(
        project_path(config["paths"]["input_metadata"]),
        partial_sequences
    )
    assert(partial_gii_sequences)

    partial_unique_sequences: set[str] = ct.get_unique_sequences(
        PARTIAL_FASTA_PATH
    )
    assert(partial_unique_sequences)

    partial_non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        PARTIAL_FASTA_PATH,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(partial_non_ambiguous_sequences)

    partial_rdrp_filtered_sequences: set[str] = partial_sequences \
        & partials_rdrp & partial_gii_sequences & partial_unique_sequences \
        & partial_non_ambiguous_sequences
    assert(partial_rdrp_filtered_sequences)

    partial_vp1_filtered_sequences: set[str] = partial_sequences \
        & partials_vp1 & partial_gii_sequences & partial_unique_sequences \
        & partial_non_ambiguous_sequences
    assert(partial_vp1_filtered_sequences)

    partial_junction_filtered_sequences: set[str] = partial_sequences \
        & partials_junction & partial_gii_sequences & partial_unique_sequences \
        & partial_non_ambiguous_sequences
    assert(partial_junction_filtered_sequences)

    partial_all_filtered_sequences = partial_rdrp_filtered_sequences \
        | partial_vp1_filtered_sequences | partial_junction_filtered_sequences

    PARTIAL_FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        PARTIAL_FASTA_PATH, 
        partial_all_filtered_sequences,
        PARTIAL_FILTERED_FASTA_PATH
    )

    # ct.typing_tool_intialise(
    #     PARTIAL_FILTERED_FASTA,
    #     "partial",
    #     batch_size=2000
    # )

#|===| END OF PARTIAL BLOCK |=================================================|
#|                                                                            |
#|===| START OF COMPLETE BLOCK |==============================================|
COMPLETE_CDS_FASTA = project_path(config["paths"]["complete_cds_fasta"])

COMPLETE_TYPING_CSV: Path | None = ct.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "complete"
)
if COMPLETE_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

COMPLETE_TYPING_CSV = project_path(COMPLETE_TYPING_CSV)

COMPLETE_GENOMIC_SPECIFICATION: Path | None = ct.get_genomic_info(
    project_path(config["paths"]["complete_full_specification"]),
    typing_information=True,
    TYPING_FILE=COMPLETE_TYPING_CSV,
    region_information=True,
    CDS_FASTA=COMPLETE_CDS_FASTA,
)

ct.extract_and_save_to_fasta(
    COMPLETE_FASTA_PATH,
    complete_filtered_sequences,
    FINAL_COMPLETE_FASTA
)
#|===| END OF COMPLETE BLOCK |================================================|
#|                                                                            |
#|===| START OF PARTIAL BLOCK |===============================================|
PARTIAL_CDS_FASTA = project_path(config["paths"]["partial_cds_fasta"])
PARTIAL_TYPING_CSV: Path | None = ct.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "partial"
)
if PARTIAL_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

PARTIAL_TYPING_CSV = project_path(PARTIAL_TYPING_CSV)

PARTIAL_SPECIFICATION: Path | None = ct.get_genomic_info(
    project_path(config["paths"]["partial_full_specification"]),
    typing_information=True,
    TYPING_FILE=PARTIAL_TYPING_CSV,
    region_information=True,
    CDS_FASTA=PARTIAL_CDS_FASTA,
)

ct.extract_and_save_to_fasta(
    PARTIAL_FASTA_PATH,
    partial_all_filtered_sequences,
    FINAL_PARTIAL_FASTA
)
#|===| END OF PARTIAL BLOCK |=================================================|

# ct.clean_up_temporary_files(TYPING_CSV)

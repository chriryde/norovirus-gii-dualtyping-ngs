import shutil
import yaml

import curation_tools as ct

from pathlib import Path

CURATE_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (CURATE_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

TEMP_TYPING_DIR = PROJECT_ROOT / Path(config["paths"]["temp_typing_dir"])
if not TEMP_TYPING_DIR.is_dir():
    print("shouldn't run this")
    RAW_FASTA = PROJECT_ROOT / Path(config["paths"]["input_fasta"])

    all_accessions: set[str] = ct.get_accessions(
        RAW_FASTA
    )

    GII_FASTA = PROJECT_ROOT / Path(config["paths"]["gii_fasta"])

    gii_sequences: set[str] = ct.get_GII_sequences(
        PROJECT_ROOT / Path(config["paths"]["input_metadata"])
    )
    assert(gii_sequences)

    complete_sequences, partial_sequences = ct.get_sequence_completeness(
        PROJECT_ROOT / Path(config["paths"]["input_metadata"])
    )
    assert(complete_sequences)

    length_filtered_sequences: set[str] = ct.get_length_filtered_sequences(
        RAW_FASTA,
        config["criteria"]["length_threshold"]
    )
    assert(length_filtered_sequences)

    unique_sequences: set[str] = ct.get_unique_sequences(
        RAW_FASTA
    )
    assert(unique_sequences)

    annotated_sequences: set[str] = ct.get_annotated_sequences(
        PROJECT_ROOT / Path(config["paths"]["input_metadata"])
    )
    assert(annotated_sequences)

    non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        RAW_FASTA,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(non_ambiguous_sequences)

    filtered_sequences: set[str] = complete_sequences & length_filtered_sequences \
        & unique_sequences & annotated_sequences & non_ambiguous_sequences & gii_sequences
    assert(filtered_sequences)

    FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        RAW_FASTA, 
        filtered_sequences,
        Path(config["paths"]["fasta_dir"]) / "filtered_sequences.fna"
    )

    FILTERED_FASTA = PROJECT_ROOT / Path(FILTERED_FASTA)

    excluded_accessions = all_accessions - filtered_sequences
    ct.export_excluded_sequences(
        Path(config["paths"]["excluded_csv"]),
        excluded_accessions=excluded_accessions,
        gii_sequences=gii_sequences,
        complete_sequences=complete_sequences,
        length_filtered_sequences=length_filtered_sequences,
        unique_sequences=unique_sequences,
        annotated_sequences=annotated_sequences,
        non_ambiguous_sequences=non_ambiguous_sequences
    )

    ct.typing_tool_intialise(
        FILTERED_FASTA
    )

    FILTERED_CDS_FASTA: Path = ct.extract_and_save_to_fasta(
        PROJECT_ROOT / Path(config["paths"]["input_cds"]), 
        filtered_sequences,
        Path(config["paths"]["fasta_dir"]) / "filtered_cds.fna"
    )

FILTERED_CDS_FASTA = PROJECT_ROOT / Path(config["paths"]["fasta_dir"]) / "filtered_cds.fna"

TYPING_CSV: Path | None = ct.typing_tool_get_results(TEMP_TYPING_DIR)
if TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    exit

src = Path(TYPING_CSV)
dst = PROJECT_ROOT / Path(config["paths"]["temp_typing_dir"]) / src.name

TYPING_CSV = Path(shutil.move(src, dst))

# ct.clean_up_temporary_files(TYPING_CSV)

GENOMIC_SPECIFICATION: Path | None = ct.get_genomic_region_info(
    TYPING_CSV,
    option="complete",
    CDS_FASTA=FILTERED_CDS_FASTA,
)

src = Path(GENOMIC_SPECIFICATION)
dst = PROJECT_ROOT / Path(config["paths"]["metadata_dir"]) / src.name

GENOMIC_SPECIFICATION = Path(shutil.move(src, dst))
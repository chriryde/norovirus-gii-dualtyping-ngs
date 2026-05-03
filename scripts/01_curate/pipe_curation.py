import pandas as pd
import typing
import subprocess
import json
import yaml

import curation_tools as ct

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

CURATE_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (CURATE_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

TEMP_TYPING_DIR = PROJECT_ROOT / Path(config["paths"]["temp_typing_dir"])
if not TEMP_TYPING_DIR.is_dir():
    print("shouldn't run this")
    GII_FASTA = PROJECT_ROOT / Path(config["paths"]["gii_fasta"])

    gii_sequences: set[str] = ct.get_GII_sequences(
        PROJECT_ROOT / Path(config["paths"]["input_fasta"]),
        PROJECT_ROOT / Path(config["paths"]["input_metadata"]),
        GII_FASTA
    )
    assert(gii_sequences)

    complete_sequences, partial_sequences = ct.get_sequence_completeness(
        PROJECT_ROOT / Path(config["paths"]["input_metadata"])
    )
    assert(complete_sequences)

    length_filtered_sequences: set[str] = ct.get_length_filtered_sequences(
        GII_FASTA,
        config["criteria"]["length_threshold"]
    )
    assert(length_filtered_sequences)

    unique_sequences: set[str] = ct.get_unique_sequences(
        GII_FASTA
    )
    assert(unique_sequences)

    annotated_sequences: set[str] = ct.get_annotated_sequences(
        PROJECT_ROOT / Path(config["paths"]["input_metadata"])
    )
    assert(annotated_sequences)

    non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        GII_FASTA,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(non_ambiguous_sequences)

    filtered_sequences: set[str] = complete_sequences & length_filtered_sequences \
        & unique_sequences & annotated_sequences & non_ambiguous_sequences
    assert(filtered_sequences)

    FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        GII_FASTA, 
        filtered_sequences,
        "filtered_sequences"
    )

    FILTERED_FASTA = PROJECT_ROOT / Path(FILTERED_FASTA)

    ct.typing_tool_intialise(
        FILTERED_FASTA
    )

    FILTERED_CDS_FASTA: Path = ct.extract_and_save_to_fasta(
        PROJECT_ROOT / Path(config["paths"]["input_cds"]), 
        filtered_sequences,
        "filtered_cds"
    )

    FILTERED_CDS_FASTA = PROJECT_ROOT / Path(FILTERED_CDS_FASTA)

TYPING_CSV: Path | None = ct.typing_tool_get_results(TEMP_TYPING_DIR)
if TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    exit

TYPING_CSV = PROJECT_ROOT / Path(TYPING_CSV)

# ct.clean_up_temporary_files(TYPING_CSV)

GENOMIC_INFORMATION: Path | None = ct.get_genomic_region_info(
    TYPING_CSV,
    option="complete",
    CDS_FASTA=FILTERED_CDS_FASTA,
)
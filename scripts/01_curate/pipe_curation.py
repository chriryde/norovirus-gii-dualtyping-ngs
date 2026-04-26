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

GII_FASTA, gii_sequences: tuple[Path, set[str]] = ct.get_GII_sequences(
    config["paths"]["input_fasta"],
    config["paths"]["input_metadata"]
)

complete_sequences, partial_sequences: tuple[set[str], set[str]] = ct.get_sequence_completeness(
    GII_FASTA
)

length_filtered_sequences: set = ct.get_length_filtered_sequences(
    GII_FASTA,
    config["parameters"]["length_threshold"]
)

unique_sequences: set = ct.get_unique_sequences(
    GII_FASTA
)

annotated_sequences: set = ct.get_annotated_sequences(
    GII_FASTA
)


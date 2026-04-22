import pandas as pd
import subprocess
import json
import csv
import yaml
from typing_crawler import run_norovirus_typing_tool

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

CURATE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (CURATE_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]

INPUT_FASTA = PROJECT_ROOT / paths_config["output_fasta"]
INPUT_METADATA = PROJECT_ROOT / paths_config["input_metadata"]

QC_DIR = PROJECT_ROOT / paths_config["qc_dir"]
METADATA_DIR = PROJECT_ROOT / paths_config["metadata_dir"]

TYPING_CSV = run_norovirus_typing_tool(INPUT_FASTA, METADATA_DIR)

def make_row():
    return {
        "BLAST_score": None,
        "length": None,
        "genotype": None,
        "p_type": None,
        "rdrp_start": None,
        "rdrp_end": None,
        "vp1_start": None,
        "vp1_end": None
    }

region_dict = defaultdict(make_row)

with open(str(TYPING_CSV, "r", encoding="utf-8")) as csvfile:
    datareader = csv.reader(csvfile)

    for row in datareader:
        accession = row["name"].split(" | ")[0]
        region_dict[accession]["BLAST"] = row["BLAST"]
        region_dict[accession]["genotype"] = row["capsid type"]
        region_dict[accession]["p_type"] = row["polymerase type"]

with open(INPUT_METADATA, "r", encoding="utf-8") as file:
    for line in file:
        record = json.loads(line)

        accession = record["accession"]
        all_accessions.add(accession)

        if record.get("isAnnotated") is not True:
            excluded_dict[accession].append("not_annotated")

        if record.get("virus", {}).get("taxId") != 122929 \
            or record.get("virus", {}).get("organismName") != "Norovirus GII":
            
            excluded_dict[accession].append("non_GII_genotype")
print(f"Process finished without errors\n")

     
import pandas as pd
import subprocess
import json
import yaml

from pathlib import Path
from Bio import SeqIO

def get_excluded(criteria, FILTERED_PATH, UNFILTERED_PATH):
    PARENT_DIR = FILTERED_PATH.parent

    UNFILTERED_TXT = PARENT_DIR / f"unfiltered_{criteria}.txt"
    FILTERED_TXT = PARENT_DIR / f"filtered_{criteria}.txt"
    EXCLUDED_TXT = PARENT_DIR / f"text_{criteria}_excluded.txt"
    EXCLUDED_FASTA = PARENT_DIR / f"temp_{criteria}_excluded.fasta"

    with open(UNFILTERED_TXT, "w") as file:
        subprocess.run(
            ["seqkit", "seq", "-n", str(UNFILTERED_PATH)],
            stdout=file,
            check=True
        )
    with open(FILTERED_TXT, "w") as file:
        subprocess.run(
            ["seqkit", "seq", "-n", str(FILTERED_PATH)]
            stdout=file,
            check=True
        )
    with open(EXCLUDED_TXT, "w") as file:
        subprocess.run(
            ["bash", "-lc", f"comm -23 <(sort {UNFILTERED_TXT}) <(sort {FILTERED_TXT})"],
            stdout=file,
            check=True
        )
    with open(EXCLUDED_FASTA, "w") as file:
        subprocess.run(
            ["seqkit", "grep", "-f", str(EXCLUDED_TXT), str(UNFILTERED_PATH)],
            stdout=file,
            check=True
        )

CURATE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (CURATE_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "curate.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]

extracted_dir = PROJECT_ROOT / paths_config["extracted_dir"]
input_fasta = PROJECT_ROOT / paths_config["input_fasta"]
input_metadata = PROJECT_ROOT / paths_config["input_metadata"]
data_dir = PROJECT_ROOT / paths_config["data_dir"]

criteria_config = config["criteria"]

filtered_length = data_dir / "filtered_length.fa"

rmlen_cmd = [
    "seqkit",
    "-m", criteria_config["length"],
    "-o", str(filtered_length),
    str(input_fasta)
]

rmdup_cmd = [
    "seqkit",
    "rmdup",
    "-s", 
    "--ignore-case"
]







excluded_accessions = []
all_accessions = []

with open(str(input_metadata), "r", encoding="utf-8") as file:
    accession = "NULL"
    row = 0
    for key, value in json.load(file).items():
        match key:
            case "accession":
                accession = key
                all_accessions.append({
                    "accession": accession
                })
            case "isAnnotated":
                if value == "False":
                    excluded_accessions.append({
                        "accession": accession,
                        "reason": "not_annotated"
                    })
            case "organismName":
                if value != "Norovirus GII":
                    excluded_accessions.append({
                        "accession": accession,
                        "reason": "non_GII_genotype"
                    })
        ind += 1
            

        
        
        

        if key == "isAnnotated" and value == "False":
            




with open(str(input_fasta)) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        accession = record.id

        all_accessions.append(accession)

        if len(record.seq) < criteria_config["length"]:
            excluded_accessions.append({
                "accession": accession,
                "reason": "too_short"
            })
        
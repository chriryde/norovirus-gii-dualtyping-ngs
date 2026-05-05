import pandas as pd
import subprocess
import json
import csv
import yaml
from typing_crawler import run_norovirus_typing_tool

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def make_row():
    return {
        "length": None,
        "BLAST_score": None,
        "genotype": None,
        "p_type": None,
        "rdrp_start": None,
        "rdrp_end": None,
        "vp1_start": None,
        "vp1_end": None
    }

CURATE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (CURATE_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]

INPUT_FASTA = PROJECT_ROOT / paths_config["output_fasta"]
INPUT_METADATA = PROJECT_ROOT / paths_config["input_annotation"]
INPUT_CDS = PROJECT_ROOT / paths_config["input_cds"]

QC_DIR = PROJECT_ROOT / paths_config["qc_dir"]
METADATA_DIR = PROJECT_ROOT / paths_config["metadata_dir"]

TYPING_CSV = run_norovirus_typing_tool(INPUT_FASTA, METADATA_DIR)

vp1_names = {
    "VP1",
    "capsid protein VP1",
    "major capsid protein",
    "major viral capsid protein",
    "major capsid protein VP1"
}

rdrp_names = {
    "RdRp",
    "RNA-dependent RNA polymerase"
}

polyprotein_names = {
     "nonstructural polyprotein"
}

typing_dict = defaultdict(make_row)

with open(TYPING_CSV, "r", encoding="utf-8") as csvfile:
    datareader = csv.DictReader(csvfile)

    for row in datareader:
        accession = row["name"] #.split(" | ")[0]
        typing_dict[accession]["length"] = row["length"]
        typing_dict[accession]["BLAST_score"] = row["BLAST score"]
        typing_dict[accession]["genotype"] = row["capsid type"]
        typing_dict[accession]["p_type"] = row["polymerase type"]

# with open(str(INPUT_FASTA)) as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         accession = record.id

for record in SeqIO.parse(INPUT_CDS, "fasta"):
        accession, coords = record.id.split(":")
        begin, end = coords.split("-")

        coding_region = record.description.replace(record.id, "", 1).strip()
        coding_region = coding_region.split("[", 1)[0].strip()

        if coding_region in vp1_names:
            typing_dict[accession]["vp1_start"] = int(begin)
            typing_dict[accession]["vp1_end"] = int(end)
            continue
            
        if coding_region in rdrp_names:
            typing_dict[accession]["rdrp_start"] = int(begin)
            typing_dict[accession]["rdrp_end"] = int(end)
            continue

        # if coding_region in polyprotein_names


# with open(INPUT_METADATA, "r", encoding="utf-8") as file:
#     for line in file:
#         record = json.loads(line)
#         accession = record["accession"]

#         for gene in record["genes"]:
#             for cds in gene["cds"]:
#                 if cds["name"] == "VP1":
#                     typing_dict[accession]["vp1_start"] = int(cds["nucleotide"]["range"][0]["begin"])
#                     typing_dict[accession]["vp1_end"] = int(cds["nucleotide"]["range"][0]["end"])
                
#                 if cds.get("name") == "RdRp":
#                     typing_dict[accession]["rdrp_start"] = int(cds["nucleotide"]["range"][0]["begin"])
#                     typing_dict[accession]["rdrp_end"] = int(cds["nucleotide"]["range"][0]["end"])

#                 for peptide in cds.get("maturePeptide", []):
#                     if peptide.get("name") == "RdRp":
#                         typing_dict[accession]["rdrp_start"] = int(peptide["nucleotide"]["range"][0]["begin"])
#                         typing_dict[accession]["rdrp_end"] = int(peptide["nucleotide"]["range"][0]["end"])

typing_df = pd.DataFrame.from_dict(typing_dict, orient="index")
typing_df.index.name = "accession"
typing_df = typing_df.reset_index()

TYPING_TSV = QC_DIR / "typing.tsv"
typing_df.to_csv(str(TYPING_TSV), sep="\t", index=False)

print(f"Process finished without errors\n")
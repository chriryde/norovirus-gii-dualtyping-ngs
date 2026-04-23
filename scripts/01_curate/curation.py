import pandas as pd
import subprocess
import json
import yaml

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

CURATE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (CURATE_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]
EXTRACTED_DIR = PROJECT_ROOT / paths_config["extracted_dir"]
INPUT_METADATA = PROJECT_ROOT / paths_config["input_metadata"]
INPUT_FASTA = PROJECT_ROOT / paths_config["input_fasta"]

QC_DIR = PROJECT_ROOT / paths_config["qc_dir"]
OUTPUT_EXCLUDED = PROJECT_ROOT / paths_config["output_excluded"]
OUTPUT_FASTA = PROJECT_ROOT / paths_config["output_fasta"]

criteria_config = config["criteria"]
excluded_dict = defaultdict(list)

print(f"Removing sequences that does not meet length criteria threshold: {criteria_config["length_threshold"]}")
try:
    rmlen_cmd = ["seqkit", "seq", "-m", str(criteria_config["length_threshold"]), "-n", str(INPUT_FASTA)]
    length_filtered = subprocess.run(
        rmlen_cmd,
        capture_output=True,
        text=True,
        check=True
    )
except subprocess.CalledProcessError as err:
    print("Command failed:")
    print("cmd", err.cmd)
    print("returncode", err.returncode)
    print("stdout", err.stdout)
    print("stderr", err.stderr)
    raise
print(f"Process finished without errors\n")

above_length_threshold = {line.split()[0] for line in length_filtered.stdout.splitlines() if line.strip()}

print(f"Removing duplicated sequences: {criteria_config["length_threshold"]}")
try:
    dededup_cmd = ["seqkit", "rmdup", "--ignore-case", "-n", str(INPUT_FASTA)]
    deduplicated_filtered = subprocess.run(
        dededup_cmd,
        capture_output=True,
        text=True,
        check=True
    )
except subprocess.CalledProcessError as err:
    print("Command failed:")
    print("cmd", err.cmd)
    print("returncode", err.returncode)
    print("stdout", err.stdout)
    print("stderr", err.stderr)
    raise
print(f"Process finished without errors\n")

not_duplicated = {line.split()[0] for line in length_filtered.stdout.splitlines() if line.strip()}

all_accessions = set()

print(f"Removing poorly annotated or non GII genotypes")
# https://stackoverflow.com/questions/50475635/loading-jsonl-file-as-json-objects
# File is JSONL not JSON
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

print(f"Removing ambiguous bases")
allowed = set("ACGT")
with open(str(INPUT_FASTA)) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        accession = record.id
        all_accessions.add(accession)

        seq = str(record.seq).upper()

        ambiguous_count = sum(1 for base in seq if base not in allowed)
        ambiguous_fraction = ambiguous_count / len(seq)
        
        if ambiguous_fraction > criteria_config["ambigious_threshold"]:
            excluded_dict[accession].append("above_ambigous_threshold")
print(f"Process finished without errors\n")

below_length_threshold = all_accessions - above_length_threshold
duplicated = all_accessions - not_duplicated

for accession in below_length_threshold:
    excluded_dict[accession].append("below_length_threshold")

for accession in duplicated:
    excluded_dict[accession].append("duplicated")

# CLUSTER_FILE = DATA_DIR / "temp_cluster_file.clstr"
# cd_hit_cmd = [
#     "cd-hit-est", "-i", str(INPUT_FASTA),
#     "-c", str(criteria_config["identity_threshold"]),
#     "-n", str(config["cd-hit"]["word_length"]),
#     "-M", str(config["parameters"]["memory"]),
#     "-T", str(config["parameters"]["threads"]),
#     "-o", str(CLUSTER_FILE)
# ]

# print(f"Starting clustering program")
# try:
#     subprocess.run(
#         cd_hit_cmd,
#         check=True
#     )
# except subprocess.CalledProcessError as err:
#     print("Command failed:")
#     print("cmd", err.cmd)
#     print("returncode", err.returncode)
#     print("stdout", err.stdout)
#     print("stderr", err.stderr)
#     raise
# print(f"Process finished without errors\n")

excluded_df = pd.DataFrame(
    columns=[
        "accession", 
        "below_length_threshold", 
        "duplicated",
        "not_annotated",
        "non_GII_genotype",
        "above_ambigous_threshold"
    ]
)

# for accession, reasons in excluded_dict.items():
#     if accession not in excluded_df.index:
#         excluded_df.loc[accession] = False

#     for reason in reasons:
#         excluded_df.loc[accession, reason] = True

# EXCLUDED_TSV = QC_DIR / "excluded.tsv"
# excluded_df.to_csv(str(EXCLUDED_TSV), sep="\t", index=False)

excluded_accessions = set(excluded_dict.keys())
kept_accessions = all_accessions - excluded_accessions

with open(INPUT_FASTA) as in_handle, open(OUTPUT_FASTA, "w") as out_handle:
    for record in SeqIO.parse(in_handle, "fasta"):
        accession = record.id

        if accession in kept_accessions:
            SeqIO.write(record, out_handle, "fasta")

print(f"Total amount of sequences: {len(all_accessions)}")
print(f"Total amount of kept accessions: {len(kept_accessions)}")
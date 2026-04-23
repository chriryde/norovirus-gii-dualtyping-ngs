import pandas as pd
import subprocess
import json
import yaml

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


ALIGN_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (ALIGN_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "alignment.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]
#EXTRACTED_DIR = PROJECT_ROOT / paths_config["extracted_dir"]
#INPUT_METADATA = PROJECT_ROOT / paths_config["input_metadata"]
INPUT_FASTA = PROJECT_ROOT / paths_config["input_fasta"]

#QC_DIR = PROJECT_ROOT / paths_config["qc_dir"]
#OUTPUT_EXCLUDED = PROJECT_ROOT / paths_config["output_excluded"]
OUTPUT_FASTA = PROJECT_ROOT / paths_config["output_fasta"]

try:
    print("running mafft")
    with open(OUTPUT_FASTA, "w") as out_f:
        mafft_cmd = ["mafft", "--auto", str(INPUT_FASTA)]
        length_filtered = subprocess.run(
            mafft_cmd,
            stdout=out_f,
            stderr=subprocess.PIPE,
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
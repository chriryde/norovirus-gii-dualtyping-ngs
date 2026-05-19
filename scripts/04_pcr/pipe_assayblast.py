import pandas as pd
import subprocess
import yaml

from pathlib import Path

def project_path(pathlike: str | Path) -> Path:
    path = Path(pathlike)
    return path if path.is_absolute() else PROJECT_ROOT / path

PCR_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (PCR_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "blast.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

DEGEN_RPIMERS = project_path(config["paths"]["primer_candidates"])

EXPANDED_PRIMERS = project_path(config["paths"]["primer_converted"])

FASTA_INPUT = project_path(config["paths"]["fasta_input"])

OUTPUT_FILE = project_path(config["paths"]["blast_output"])

ANALYZE_OUTPUT = project_path(config["paths"]["assay_analysis"])

expand_primers = project_path(config["paths"]["expand_script"])

cmd = [
    "python",
    expand_primers,
    "-i", str(DEGEN_RPIMERS),
    "-o", str(EXPANDED_PRIMERS),
]

subprocess.run(cmd, check=True)

cmd = [
    "assay_blast",
    str(FASTA_INPUT),
    "-q", str(EXPANDED_PRIMERS),
    "-o", str(OUTPUT_FILE),
    "--mismatch", "3"
]

result = subprocess.run(
    cmd,
    text=True,
    capture_output=True,
)

print("STDOUT:")
print(result.stdout)

print("STDERR:")
print(result.stderr)

cmd_analyze = [
    "assay_analyze",
    str(OUTPUT_FILE),
    "-o", str(ANALYZE_OUTPUT),
    "--mismatch", "3",
    "--only-primer",
    "--distance", "8000"
]

result = subprocess.run(
    cmd_analyze,
    check=True,
    text=True,
    capture_output=True,
)

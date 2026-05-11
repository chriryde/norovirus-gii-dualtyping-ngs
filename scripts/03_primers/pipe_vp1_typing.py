from pathlib import Path
import primer_tools as pt
from Bio import SeqIO
import yaml

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "typing.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)



INPUT_FASTA = PROJECT_ROOT / Path(config["paths"]["vp1_partial"])
EXTRACTED_FASTA = PROJECT_ROOT / Path(config["paths"]["output_fasta_vp1_amplicon"])
OUTPUT_FASTA_DIR = PROJECT_ROOT / Path(config["paths"]["output_partials_fasta_dir"])

TYPING_DIR = PROJECT_ROOT / Path(config["paths"]["web_crawler_dir"])

#kommentera in detta när vi ska köra ett nytt jobb
# pt.typing_tool_intialise(
#     EXTRACTED_FASTA,
#     "primer_pair_3_partial_vp1_run"
# )

VP1_TYPING_CSV: Path | None = pt.typing_tool_get_results(
    TYPING_DIR,
    "primer_pair_3_partial_vp1_run"
)
if VP1_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit


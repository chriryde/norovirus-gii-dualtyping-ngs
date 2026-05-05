from pathlib import Path
import primer_tools as pt
from Bio import SeqIO
import yaml

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "amplicon_typing.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

#gör om sökvägen sen
# COORD_VP1_PATH = PROJECT_ROOT / "data" / "03_primer_evaluation" /"novel" / "varvamp_outputs"/"test_output"/"filtered_primers.tvs"
# COORD_RDRP_PATH = PROJECT_ROOT / "data" / "03_primer_evaluation" /"novel" / "varvamp_outputs"/"test_output"/"filtered_primers.tvs"

TEST_PRIMERS = PROJECT_ROOT / Path(config["paths"]["input_dir"]) / "test_output" / "filtered_primers.tsv"

INPUT_FASTA = PROJECT_ROOT / Path(config["paths"]["input_fasta"])
EXTRACTED_FASTA = PROJECT_ROOT / Path(config["paths"]["extracted_fasta"])
OUTPUT_FASTA_DIR = PROJECT_ROOT / Path(config["paths"]["output_fasta_dir"])

TEMP_TYPING_DIR = PROJECT_ROOT / Path(config["paths"]["temporary_web_crawler_dir"])

COMBINED_PRIMERS: Path = PROJECT_ROOT / Path(
    pt.extract_amplicons(
        TEST_PRIMERS,
        6,
        TEST_PRIMERS,
        6,
        INPUT_FASTA,
        OUTPUT_FASTA_DIR
    )
)

# pt.typing_tool_intialise(
#     COMBINED_PRIMERS,
#     "test"
# )

COMBINED_TYPING_CSV: Path | None = pt.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "test"
)
if COMBINED_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

# regions = []

# print("starting amplicon typing")

# with open(COORD_RDRP_PATH, "r", encoding="utf-8") as file:
#     for line in file:
#         if line.strip() == "":
#             continue

#         cols = line.strip().split()
#         start = int(cols[1])
#         #end = int(cols[2])
#         regions.append((start))

# with open(COORD_VP1_PATH, "r", encoding="utf-8") as file:
#     for line in file:
#         if line.strip() == "":
#             continue

#         cols = line.strip().split()
#         #start = int(cols[1])
#         end = int(cols[2])
#         regions.append((end))

# print("Regions:", regions)

# print("Extracting:", start, stop)
# with open(EXTRACTED_FASTA, "w") as out:
#     for record in SeqIO.parse(INPUT_FASTA, "fasta"):
#         extracted_seq = record.seq[start:stop]

#         out.write(f">{record.id}_region_{overall_start}_{overall_end}\n")
#         out.write(str(extracted_seq) + "\n")

# print(f"Amplicons saved in {EXTRACTED_FASTA}")
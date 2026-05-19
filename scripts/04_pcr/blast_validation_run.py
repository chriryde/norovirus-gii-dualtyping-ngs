from pathlib import Path
import yaml

import BLAST_validation as bv

PCR_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (PCR_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "blast.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]


ASSAY_DETAILS = PROJECT_ROOT / paths_config["assay_details"]

AMPLICON_CSV = PROJECT_ROOT / paths_config["amplicon_csv"]

SPECIFICATION_CSV = PROJECT_ROOT / paths_config["specification_csv"]

AMPLICON_SPECIFICATION_CSV = PROJECT_ROOT / paths_config["amplicon_specification_csv"]


genome_dict, primer_interaction_dict = bv.check_amplicons(ASSAY_DETAILS)

# print(genome_dict)
print(primer_interaction_dict)

bv.passed_amplicons(genome_dict, AMPLICON_CSV)

bv.passed_types(AMPLICON_CSV, SPECIFICATION_CSV, AMPLICON_SPECIFICATION_CSV)

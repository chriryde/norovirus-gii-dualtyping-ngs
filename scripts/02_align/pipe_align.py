import alignment_tools as at

## från global alignment:
import pandas as pd
import subprocess
import json
import yaml

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


ALIGN_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (ALIGN_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "alignment.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]
##########


## runs global alignment, stores raw msa file + global log file
# GLOBAL_FASTA, GLOBAL_LOG = at.alignment(
#     PROJECT_ROOT / paths_config["final_global"], 
#     PROJECT_ROOT / paths_config["output_global"],
#     PROJECT_ROOT /paths_config['global_log']
# ) 

# # # ## runs rdrp alignment, stores raw msa file + rdrp log file
RDRP_ALIGN, RDRP_LOG = at.alignment(
    PROJECT_ROOT /paths_config["polymerase_correct"], 
    PROJECT_ROOT /paths_config["output_msa_rdrp"],
    PROJECT_ROOT /paths_config['rdrp_log']
) 


# # # ## runs vp1 alignment, stores raw msa file + vp1 log file
VP1_ALIGN, VP1_LOG = at.alignment(
    PROJECT_ROOT /paths_config["vp1"], 
    PROJECT_ROOT /paths_config["output_msa_vp1"],
    PROJECT_ROOT /paths_config['vp1_log']
) 


# # # ## runs rdrp+vp1 alignment, stores raw msa file + rdrp+vp1 log file
# RDRP_VP1_ALIGN, RDRP_VP1_LOG = at.alignment(
#     PROJECT_ROOT / paths_config["final_rdrp_vp1"], 
#     PROJECT_ROOT / paths_config["output_rdrp_vp1"],
#     PROJECT_ROOT / paths_config['rdrp_vp1_log']
# ) 

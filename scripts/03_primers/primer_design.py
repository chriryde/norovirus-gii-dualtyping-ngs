import yaml
import subprocess
import csv
import primer_tools as pt


from collections import defaultdict
from pathlib import Path


#paths to files outside this directory, fetched from config 
PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (PRIMER_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "fetch_to_primer.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]
INPUT_FASTA = PROJECT_ROOT / paths_config["input_fasta"]
INPUT_FASTA_RDRP = PROJECT_ROOT / paths_config["rdrp_aln"]
INPUT_FASTA_VP1 = PROJECT_ROOT / paths_config["vp1_aln"]

#QC_DIR = PROJECT_ROOT / paths_config["qc_dir"]
#OUTPUT_EXCLUDED = PROJECT_ROOT / paths_config["output_excluded"]
OUTPUT_DIR = PROJECT_ROOT / paths_config["output_dir"]
OUTPUT_DIR_NAME_COMPLETE = PROJECT_ROOT / paths_config["output_dir"] / "test_output"
OUTPUT_DIR_NAME_RDRP = PROJECT_ROOT / paths_config["output_dir"] / "test_output_rdrp"
OUTPUT_DIR_NAME_VP1 = PROJECT_ROOT / paths_config["output_dir"] / "test_output_vp1"

TSV_PRIMERS = PROJECT_ROOT / paths_config["output_dir"] / "test_output" / "primers.tsv"
TSV_FILTERED_PRIMERS = PROJECT_ROOT / paths_config["output_dir"] / "test_output" / "filtered_primers.tsv"

TSV_PRIMERS_RDRP = PROJECT_ROOT / paths_config["output_dir"] / "test_output_rdrp" / "primers.tsv"
TSV_FILTERED_PRIMERS_RDRP = PROJECT_ROOT / paths_config["output_dir"] / "test_output_rdrp" / "filtered_primers.tsv"

TSV_PRIMERS_VP1 = PROJECT_ROOT / paths_config["output_dir"] / "test_output_vp1" / "primers.tsv"
TSV_FILTERED_PRIMERS_VP1 = PROJECT_ROOT / paths_config["output_dir"] / "test_output_vp1" / "filtered_primers.tsv"

REFERENCE_LIBRARY = PROJECT_ROOT / paths_config["reference_library"] / "off_targets"

REFERENCE_FASTA = PROJECT_ROOT / paths_config["final_complete_fna"]

#fetching parameters for varVamp from config
config_varvamp = config["varvamp"]
scheme = config_varvamp["scheme"]
opt_length = config_varvamp["opt_length_complete"]
max_length = config_varvamp["max_length_complete"]
opt_length_rdrp = config_varvamp["opt_length_rdrp"]
max_length_rdrp = config_varvamp["max_length_rdrp"]
#n_ambig = config_varvamp["n_ambig"]

#Run varVamp through linux with reference library

# OUTPUT_COMPLETE_DIR = pt.varvamp(scheme, opt_length, max_length,
#             REFERENCE_LIBRARY, INPUT_FASTA, OUTPUT_DIR, 'complete')

# OUTPUT_COMPLETE_DIR = pt.varvamp_fast(scheme, opt_length, max_length,
#                         INPUT_FASTA, OUTPUT_DIR, 'complete')

# PATH_TO_PRIMER_TSV = OUTPUT_COMPLETE_DIR / 'primers.tsv'
# PATH_TO_PRIMER_BED = OUTPUT_COMPLETE_DIR / 'primers.bed'
# pt.correct_primer_position(REFERENCE_FASTA, PATH_TO_PRIMER_TSV, PATH_TO_PRIMER_BED)

# TSV_PRIMERS_COMPLETE = OUTPUT_COMPLETE_DIR / "primers.tsv"
# pt.filter(TSV_PRIMERS_COMPLETE, OUTPUT_COMPLETE_DIR)


#varVamp for RDRP

OUTPUT_RDRP_DIR = pt.varvamp(scheme, opt_length_rdrp, max_length_rdrp, REFERENCE_LIBRARY,
            INPUT_FASTA_RDRP, OUTPUT_DIR, 'rdrp_initial')

PATH_TO_PRIMER_TSV = OUTPUT_RDRP_DIR / 'primers.tsv'
PATH_TO_PRIMER_BED = OUTPUT_RDRP_DIR / 'primers.bed'
pt.correct_primer_position(REFERENCE_FASTA, PATH_TO_PRIMER_TSV, PATH_TO_PRIMER_BED)

# TSV_PRIMERS_RDRP = OUTPUT_RDRP_DIR / "primers.tsv"
# pt.filter(TSV_PRIMERS_RDRP, OUTPUT_RDRP_DIR)

#varVamp for VP1

#Run varVamp through linux with reference library
OUTPUT_VP1_DIR = pt.varvamp(scheme, opt_length_rdrp, max_length_rdrp,
            REFERENCE_LIBRARY, INPUT_FASTA_VP1, OUTPUT_DIR, 'vp1_initial')

PATH_TO_PRIMER_TSV = OUTPUT_VP1_DIR / 'primers.tsv'
PATH_TO_PRIMER_BED = OUTPUT_VP1_DIR / 'primers.bed'
pt.correct_primer_position(REFERENCE_FASTA, PATH_TO_PRIMER_TSV, PATH_TO_PRIMER_BED)

# TSV_PRIMERS_VP1 = OUTPUT_VP1_DIR / "primers.tsv"
# pt.filter(TSV_PRIMERS_VP1, OUTPUT_VP1_DIR)


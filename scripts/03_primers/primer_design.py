import yaml
import subprocess
import json
import csv

from collections import defaultdict
from pathlib import Path


PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (PRIMER_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "fetch_to_primer.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]
INPUT_FASTA = PROJECT_ROOT / paths_config["input_fasta"]

#QC_DIR = PROJECT_ROOT / paths_config["qc_dir"]
#OUTPUT_EXCLUDED = PROJECT_ROOT / paths_config["output_excluded"]
OUTPUT_DIR = PROJECT_ROOT / paths_config["output_dir"]
OUTPUT_DIR_NAME = PROJECT_ROOT / paths_config["output_dir"] / "test_output"
TSV_PRIMERS = PROJECT_ROOT / paths_config["output_dir"] / "test_output" / "primers.tsv"
TSV_FILTERED_PRIMERS = PROJECT_ROOT / paths_config["output_dir"] / "test_output" / "filtered_primers.tsv"
REFERENCE_LIBRARY = PROJECT_ROOT / paths_config["reference_library"] / "off_targets"

config_varvamp = config["varvamp"]
scheme = config_varvamp["scheme"]
opt_length = config_varvamp["opt_length"]
max_length = config_varvamp["max_length"]
#n_ambig = config_varvamp["n_ambig"]


print('startar primerdesign')
try:
    cmd = ['varvamp', str(scheme), '-ol', str(opt_length), '-ml', str(max_length),
            #'-db', str(REFERENCE_LIBRARY), 
            str(INPUT_FASTA), str(OUTPUT_DIR_NAME)]
    deduplicated_filtered = subprocess.run(
        cmd,
        check=True
    )
except subprocess.CalledProcessError as err:
    print("Command failed:")
    print("cmd", err.cmd)
    print("returncode", err.returncode)
    print("stdout", err.stdout)
    print("stderr", err.stderr)
    raise
print('avslutar primerdesign') 


print('removing primers outside bp 4000-6500')

all_primers = set()
excluded_dict = set()
counter = 0

with open(TSV_PRIMERS, 'r') as file:
    tsv_reader = csv.reader(file, delimiter='\t')
    for row in tsv_reader:
        if counter != 0:
            all_primers.add(tuple(row))
            start = int(row[5])
            if type(start) is int:
                if start < 4500 or start > 6000:
                    excluded_dict.add(tuple(row))
        counter += 1
    
    curated_primers = all_primers - excluded_dict
    
    header = ('amlicon_name', 'amplicon_length', 'primer_name', 'primer_name_all_primers', 
               'pool', 'start', 'stop', 'seq', 'size', 'gc_best', 'temp_best', 'mean_gc', 
                 'mean_temp', 'penalty', 'off_target_amplicons')
    
    with open(TSV_FILTERED_PRIMERS, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(header)
        writer.writerows(sorted(curated_primers))
    print(f'Total amount of potential primers after filtration: {len(curated_primers)}')    

print('procedure proceded without errors')
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
            '-db', str(REFERENCE_LIBRARY), str(INPUT_FASTA), str(OUTPUT_DIR_NAME)]
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
# def load_config(path="fetch_to_primer.yml"):
#     with open(path) as f:
#         return yaml.safe_load(f)
# x

# #name är det som vi skriver i vår config fil tex global_aln, skapar en egen mapp för varje namn 
# def run_varvamp(name, msa_path, output_dir, varvamp_config):
#     run_out = output_dir / name
#     run_out.mkdir(parents=True, exist_ok=True)

#     scheme = varvamp_config["scheme"]
#     cmd = ["varvamp", scheme, str(msa_path), str(run_out)]

#     # Optional overrides — only added if present and not commented out in config, ändrar våra parametrar i varvamp
#     for param in ["opt-length", "max-length", "threshold", "n-ambig", "database"]:
#         value = varvamp_config.get(param)
#         if value is not None:
#             cmd += [f"--{param}", str(value)]

#     print(f"\n[varvamp] Running on {name}...")
#     print("  CMD:", " ".join(cmd))

#     result = subprocess.run(cmd, capture_output=True, text=True)

#     if result.returncode != 0:
#         print(f"[ERROR] varVAMP failed for {name}:\n{result.stderr}")
#     else:
#         print(f"[done]  Results in {run_out}/")


# def main():
#     config = load_config()

#     output_dir = Path(config["output"]["dir"])
#     varvamp_config = config["varvamp"]
#     alignments = config["alignments"]  # → {"global_aln": "...", "local_aln": "...", "2local_aln": "..."}

#     for name, path in alignments.items():
#         msa_path = Path(path)

#         if not msa_path.exists():
#             print(f"[SKIP] {name}: file not found at {msa_path}")
#             continue

#         run_varvamp(name, msa_path, output_dir, varvamp_config)


# if __name__ == "__main__":
#     main()
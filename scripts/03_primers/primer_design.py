import yaml
import subprocess
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
TVS_PRIMERS = PROJECT_ROOT / paths_config["output_dir"] / "test_output" / 'primers.tvs'

config_varvamp = config["varvamp"]
scheme = config_varvamp["scheme"]
opt_length = config_varvamp["opt_length"]
max_length = config_varvamp["max_length"]
#n_ambig = config_varvamp["n_ambig"]


print('startar primerdesign')
try:
    cmd = ['varvamp', str(scheme), '-ol', str(opt_length), '-ml', str(max_length), str(INPUT_FASTA), str(OUTPUT_DIR_NAME)]
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

with open(TVS_PRIMERS, 'r') as file:
    for line in file:

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
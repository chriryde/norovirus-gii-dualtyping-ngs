import yaml
import subprocess
from pathlib import Path

def load_config(path="fetch_to_primer.yml"):
    with open(path) as f:
        return yaml.safe_load(f)

def run_varvamp(name, msa_path, output_dir, varvamp_config):
    run_out = output_dir / name
    run_out.mkdir(parents=True, exist_ok=True)

    scheme = varvamp_config["scheme"]
    cmd = ["varvamp", scheme, str(msa_path), str(run_out)]

    # Optional overrides — only added if present and not commented out in config
    for param in ["opt-length", "max-length", "threshold", "n-ambig", "database"]:
        value = varvamp_config.get(param)
        if value is not None:
            cmd += [f"--{param}", str(value)]

    print(f"\n[varvamp] Running on {name}...")
    print("  CMD:", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"[ERROR] varVAMP failed for {name}:\n{result.stderr}")
    else:
        print(f"[done]  Results in {run_out}/")


def main():
    config = load_config()

    output_dir = Path(config["output"]["dir"])
    varvamp_config = config["varvamp"]
    alignments = config["alignments"]  # → {"global_aln": "...", "local_aln": "...", "2local_aln": "..."}

    for name, path in alignments.items():
        msa_path = Path(path)

        if not msa_path.exists():
            print(f"[SKIP] {name}: file not found at {msa_path}")
            continue

        run_varvamp(name, msa_path, output_dir, varvamp_config)


if __name__ == "__main__":
    main()
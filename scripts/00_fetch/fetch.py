import subprocess
import yaml

from zipfile import ZipFile
from pathlib import Path

# Get the parent directory of the file
# https://stackoverflow.com/questions/3430372/how-do-i-get-the-full-path-of-the-current-files-directory
FETCH_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (FETCH_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "fetch.yml"


# https://www.geeksforgeeks.org/python/reading-and-writing-yaml-file-in-python/
with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

# Base for query command
cmd = [
    "datasets", 
    "download", 
    "virus", 
    "genome", 
    "taxon", config["taxon"]
] # Extract configuration for NCBIdataset virus genome query

download_config = config["download"]
paths_config = config["paths"]

if download_config.get("host"):
    cmd.extend(["--host", download_config["host"]])

if download_config.get("complete_only", False):
    cmd.append("--complete-only")

include_config = download_config["include"]

for include_name, enabled in include_config.items():
    if enabled:
        cmd.extend(["--include", include_name])

download_dir = PROJECT_ROOT / paths_config["download_dir"]
download_zip = PROJECT_ROOT / paths_config["download_zip"]
extracted_dir = PROJECT_ROOT / paths_config["extracted_dir"]

download_dir.mkdir(parents=True, exist_ok=True)
extracted_dir.mkdir(parents=True, exist_ok=True)

cmd.extend(["--filename", str(download_zip)])

print("Fetching from NCBIdatasets:")
print(" ".join(cmd))

subprocess.run(cmd, check=True)

if not download_zip.exists():
    raise FileNotFoundError(f"\nThe requested zip could not be downloaded: {download_zip}\n")

print(f"Extracting {download_zip} to {extracted_dir}\n")

# https://www.geeksforgeeks.org/python/unzipping-files-in-python/
with ZipFile(str(download_zip), 'r') as zObject:
    zObject.extractall(path=str(extracted_dir))
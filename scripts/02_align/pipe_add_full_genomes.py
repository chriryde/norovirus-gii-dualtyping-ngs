import alignment_tools as at
import cd_hit as cd
import partial as pa
from pathlib import Path
import yaml
import csv
import include_more_sequences_global as imsg

ALIGN_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (ALIGN_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "alignment.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]


# Get accessions
full_genomes_from_partial = imsg.global_accessions(
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["partial_metadata"]
)

# Retrieve sequences
cd.get_sequences(
    PROJECT_ROOT / paths_config["input_partial"],
    PROJECT_ROOT / paths_config["added_full_fasta"],
    full_genomes_from_partial
)

FULL_GENOMES_FROM_PARTIAL_TSV = PROJECT_ROOT / paths_config["full_genomes_from_partial_tsv"]

with open(FULL_GENOMES_FROM_PARTIAL_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['accession'])
        for accession in full_genomes_from_partial:
            writer.writerow([accession])


# Merge fasta files
pa.merge_partial_regions(
    PROJECT_ROOT / paths_config["final_global"],
    PROJECT_ROOT / paths_config["added_full_fasta"],
    PROJECT_ROOT / paths_config["final_global_with_added"]
)
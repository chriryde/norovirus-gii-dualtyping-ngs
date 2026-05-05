import csv
import yaml
from pathlib import Path

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "fetch_to_primer.yml"
with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

TEMPORARY_WEB_CRAWLER_DATA_CSV = PROJECT_ROOT / Path(config["paths"]["web_crawler_data"]) / "test_typing_data.csv"
REFERENCE_SEQUENCE_CSV = Path(config["paths"]["metadata_dir"]) 

all_typed_sequences = {}
all_reference_sequences = {}
matches = []
missmatches = []



with open(TEMPORARY_WEB_CRAWLER_DATA_CSV, 'r') as file:
    csv_reader = csv.reader(file, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        ascension_nr = str(row[0])
        polymerase_type = str(row[8])
        polymerase_subtype = str(row[14])
        #lägg till alla typningar i set
        all_typed_sequences[ascension_nr] = [polymerase_type, polymerase_subtype]


with open(REFERENCE_SEQUENCE_CSV, 'r') as file:
    csv_reader = csv.reader(file, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        ascension_nr = str(row[0])
        polymerase_type = str(row[8])
        polymerase_subtype = str(row[14])
        #lägg till alla typningar i set
        all_typed_sequences[ascension_nr] = [polymerase_type, polymerase_subtype]

for typed_sequence, reference_sequence in zip(all_typed_sequences, all_reference_sequences):
    if typed_sequence == reference_sequence:
        matches.append(typed_sequence)
    else:
        missmatches.append(typed_sequence)

print(f'Ratio of sequences correctly typed: {len(matches) / len(all_typed_sequences)}')
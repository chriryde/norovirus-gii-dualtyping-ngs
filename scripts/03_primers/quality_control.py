import csv
import yaml
from pathlib import Path

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "fetch_to_primer.yml"
with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

TEMPORARY_WEB_CRAWLER_DATA_CSV = PROJECT_ROOT / Path(config["paths"]["web_crawler_data"]) / "test_typing_data.csv"
REFERENCE_SEQUENCE_CSV = Path(config["paths"]["complete_sequence"])

all_typed_sequences = {}
all_reference_sequences = {}
matches = []
missmatches = []
matches_rdrp = []
missmatches_rdrp = []
matches_vp1= []
missmatches_vp1 = []



with open(TEMPORARY_WEB_CRAWLER_DATA_CSV, 'r') as file:
    csv_reader = csv.DictReader(file, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        ascension_nr = str(row['name']).split("_region")[0]
        polymerase_type = str(row['polymerase type'])
        polymerase_subtype = str(row['capsid type'])
        #lägg till alla typningar i set
        all_typed_sequences[ascension_nr] = [polymerase_type, polymerase_subtype]


with open(REFERENCE_SEQUENCE_CSV, 'r') as file:
    csv_reader = csv.DictReader(file, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        ascension_nr = str(row['accession'])
        polymerase_type = str(row['p_type'])
        polymerase_subtype = str(row['genotype'])
        #lägg till alla typningar i set
        all_reference_sequences[ascension_nr] = [polymerase_type, polymerase_subtype]


#check matches for both
for typed_sequence in all_typed_sequences.keys():
    if typed_sequence in all_reference_sequences:
        if all_typed_sequences[typed_sequence] == 'Could not assign':
            continue
        if all_typed_sequences[typed_sequence] == all_reference_sequences[typed_sequence]:
            matches.append(typed_sequence)
        else:
            missmatches.append(typed_sequence)
    else:
        print(f'error: accession number {typed_sequence} does not exist in reference file')
        continue

#check matches for rdrp
for typed_sequence in all_typed_sequences.keys():
    if typed_sequence in all_reference_sequences:
        if all_typed_sequences[typed_sequence] == 'Could not assign':
            continue
        if all_typed_sequences[typed_sequence][0] == all_reference_sequences[typed_sequence][0]:
            matches_rdrp.append(typed_sequence)
        else:
            missmatches_rdrp.append(typed_sequence)
    else:
        continue

#check for matches in vp1
for typed_sequence in all_typed_sequences.keys():
    if typed_sequence in all_reference_sequences:
        if all_typed_sequences[typed_sequence] == 'Could not assign':
            continue
        if all_typed_sequences[typed_sequence][1] == all_reference_sequences[typed_sequence][1]:
            matches_vp1.append(typed_sequence)
        else:
            missmatches_vp1.append(typed_sequence)
    else:
        continue

print(f'Ratio of sequences correctly typed: {len(matches) / len(all_typed_sequences)}')
print(f'Ratio of rdrp:s correctly typed: {len(matches_rdrp) / len(all_typed_sequences)}')
print(f'Ratio of vp1:s correctly typed: {len(matches_vp1) / len(all_typed_sequences)}')
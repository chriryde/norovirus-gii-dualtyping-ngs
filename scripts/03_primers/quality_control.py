import csv
import yaml
from pathlib import Path
from collections import defaultdict

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "fetch_to_primer.yml"
with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

##paths for longer amplicon
#WEB_CRAWLER_DATA_CSV_VP1 = PROJECT_ROOT / Path(config["paths"]["web_crawler_data"]) / "primer_pair_3_rdrp_run_typing_results_20260511_142616.csv"
WEB_CRAWLER_DATA_CSV = PROJECT_ROOT / Path(config["paths"]["web_crawler_data"]) / "primer_pair_3_partial_vp1_run_typing_results_20260511_144109.csv"
REFERENCE_SEQUENCE_CSV = Path(config["paths"]["merged_reference"])

##paths for short amplicon
# WEB_CRAWLER_DATA_CSV = PROJECT_ROOT / Path(config["paths"]["web_crawler_data_short"])
# REFERENCE_SEQUENCE_CSV = PROJECT_ROOT / Path(config["paths"]["complete_specification"])


# PRIMER_DIR = Path(__file__).resolve().parent
# PROJECT_ROOT = PRIMER_DIR.parent.parent

# CONFIG_PATH = PROJECT_ROOT / "config" / "fetch_to_primer.yml"
# with open(CONFIG_PATH, "r", encoding="utf-8") as f:
#     config = yaml.safe_load(f)

# WEB_CRAWLER_DATA_CSV = PROJECT_ROOT / Path(config["paths"]["web_crawler_data_short"]) 
# REFERENCE_SEQUENCE_CSV = Path(config["paths"]["complete_specification"])

all_typed_sequences = {}
all_reference_sequences = {}

matches = []
missmatches_detailed = defaultdict(list)
matches_rdrp = []
missmatches_rdrp_detailed = defaultdict(list)
matches_vp1= []
missmatches_vp1_detailed = defaultdict(list)



with open(WEB_CRAWLER_DATA_CSV, 'r') as file:
    csv_reader = csv.DictReader(file, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        ascension_nr = str(row['name']).split("_region")[0]
        polymerase_type = str(row['polymerase type']).split(" (")[0]
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

number_of_non_existant_accession_numbers = 0
#check matches for both
for typed_sequence in all_typed_sequences.keys():
    if typed_sequence in all_reference_sequences:
        if all_typed_sequences[typed_sequence] == 'Could not assign':
            continue
        
        typed = all_typed_sequences[typed_sequence]
        reference = all_reference_sequences[typed_sequence]

        if all_typed_sequences[typed_sequence] == all_reference_sequences[typed_sequence]:
            matches.append(typed_sequence)
        else:
            mismatch_key = f"Expected: {typed}  Actual: {reference}"
            missmatches_detailed[mismatch_key].append(typed_sequence)
    else:
        number_of_non_existant_accession_numbers += 1
        print(f'error: accession number {typed_sequence} does not exist in reference file')        
        continue

print(matches)

#check matches for rdrp
for typed_sequence in all_typed_sequences.keys():
    if typed_sequence in all_reference_sequences:
        if all_typed_sequences[typed_sequence] == 'Could not assign':
            continue
        
        typed_rdrp = all_typed_sequences[typed_sequence][0]
        reference_rdrp = all_reference_sequences[typed_sequence][0]
        
        if all_typed_sequences[typed_sequence][0] == all_reference_sequences[typed_sequence][0]:
            matches_rdrp.append(typed_sequence)
        else:
            mismatch_key = f"RDRP typed: {typed_rdrp}, reference type: {reference_rdrp}"
            missmatches_rdrp_detailed[mismatch_key].append(typed_sequence)
    else:
        continue

#check for matches in vp1
for typed_sequence in all_typed_sequences.keys():
    if typed_sequence in all_reference_sequences:
        if all_typed_sequences[typed_sequence] == 'Could not assign':
            continue
        
        typed_vp1 = all_typed_sequences[typed_sequence][1]
        reference_vp1 = all_reference_sequences[typed_sequence][1]
        
        if all_typed_sequences[typed_sequence][1] == all_reference_sequences[typed_sequence][1]:
            matches_vp1.append(typed_sequence)
        else:
            #print(all_typed_sequences[typed_sequence][1])
            #print(all_reference_sequences[typed_sequence][1])
            mismatch_key = f"VP1 Typed: {typed_vp1}  Reference typed: {reference_vp1}"
            missmatches_vp1_detailed[mismatch_key].append(typed_sequence)
    else:
        continue

total = len(all_typed_sequences) 

print(f"\nTotal results:")
print(f"  Correctly typed (both): {len(matches)}/{total} ({100*len(matches)/total:.1f}%)")
print(f"  RDRP correct: {len(matches_rdrp)}/{total} ({100*len(matches_rdrp)/total:.1f}%)")
print(f"  VP1 correct: {len(matches_vp1)}/{total} ({100*len(matches_vp1)/total:.1f}%)")

print(f"\n" + "=" * 70)
print("Detailed errors - RDRP (Polymerase Type)")
print("=" * 70)
if missmatches_rdrp_detailed:
    for error_type, accessions in sorted(missmatches_rdrp_detailed.items(), key=lambda x: -len(x[1])):
        print(f"\n{error_type}")
        print(f" Total occurances : {len(accessions)}")
        
else:
    print("No errors classifying RDRP")

print(f"\n" + "=" * 70)
print("Detailed errors - VP1 (Capsid Type)")
print("=" * 70)

if missmatches_vp1_detailed:
    for error_type, accessions in sorted(missmatches_vp1_detailed.items(), key=lambda x: -len(x[1])):
        print(f"\n{error_type}")
        print(f" Total occurances : {len(accessions)}")
else:
    print("No errors classifying VP1")

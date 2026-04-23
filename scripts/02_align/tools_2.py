import csv
from pathlib import Path
import yaml
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq




ALIGN_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (ALIGN_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "alignment.yml"


with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]

INPUT_METADATA = PROJECT_ROOT / paths_config["input_metadata"]
INPUT_FASTA = PROJECT_ROOT / paths_config['input_fasta']
RDRP_FASTA = PROJECT_ROOT / paths_config['rdrp_fasta']

df = pd.read_csv(INPUT_METADATA)

rdrp_dict = {}

for index, row in df.iterrows():
    rdrp_dict[row['accession']] = (row['rdrp_start'], row['rdrp_end'])

# selected_col = df.loc[:, ['accession', 'rdrp_start', 'rdrp_end']]
# print(selected_col)
# print(rdrp_dict)



# with open(INPUT_FASTA) as in_handle, open(OUTPUT_FASTA, "w") as out_handle:
#     for record in SeqIO.parse(in_handle, "fasta"):
#         accession = record.id

#         if accession in kept_accessions:
#             SeqIO.write(record, out_handle, "fasta")

with open(str(INPUT_FASTA)) as in_handle, open(RDRP_FASTA, "w") as out_handle:
    for record in SeqIO.parse(in_handle, "fasta"):
        accession = record.id
        seq = str(record.seq).upper()

        start = rdrp_dict[accession][0]
        end = rdrp_dict[accession][1]
        region = seq[start: end+1]
        # print(record)
        # print(region)
        print()

        SeqIO.write(Seq(region), out_handle, "fasta")








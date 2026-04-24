import csv
from pathlib import Path
import yaml
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




ALIGN_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = (ALIGN_DIR.parent).parent
CONFIG_PATH = PROJECT_ROOT / "config" / "alignment.yml"


with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]

INPUT_METADATA = PROJECT_ROOT / paths_config["input_metadata"]
INPUT_FASTA = PROJECT_ROOT / paths_config['input_fasta']
RDRP_FASTA = PROJECT_ROOT / paths_config['rdrp_fasta']
VP1_FASTA = PROJECT_ROOT / paths_config['vp1_fasta']

df = pd.read_csv(INPUT_METADATA)

rdrp_dict = {}
vp1_dict = {}

# check that all names match the csv file
for index, row in df.iterrows():
    rdrp_dict[row['accession']] = (row['rdrp_start'], row['rdrp_end'])
    vp1_dict[row['accession']] = (row['vp1_start'], row['vp1_end'])



with open(str(INPUT_FASTA)) as in_handle, open(RDRP_FASTA, "w") as out_handle1, open(VP1_FASTA, "w") as out_handle2:
    for record in SeqIO.parse(in_handle, "fasta"):
        accession = record.id
        seq = str(record.seq).upper()

        rdrp_start = rdrp_dict[accession][0]
        rdrp_end = rdrp_dict[accession][1]
        rdrp_region = seq[rdrp_start: rdrp_end+1]
       
        rdrp_record = SeqRecord(
            Seq(rdrp_region),
            id = accession,
            name = accession,
            description = record.description
        )

        vp1_start = vp1_dict[accession][0]
        vp1_end = vp1_dict[accession][1]
        vp1_region = seq[vp1_start : vp1_end+1]

        vp1_record = SeqRecord(
            Seq(vp1_region),
            id = accession,
            name = accession,
            description = record.description
        )
       
        SeqIO.write(rdrp_record, out_handle1, "fasta")
        SeqIO.write(vp1_record, out_handle2, "fasta")









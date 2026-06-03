import alignment_tools as at
import cd_hit as cd
import partial as pa
from pathlib import Path
import yaml
import csv
import pandas as pd

ALIGN_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (ALIGN_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "alignment.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]


### load metadata for complete and partial sequences
df_complete = pd.read_csv(str(PROJECT_ROOT / paths_config["input_metadata"]))
df_partial = pd.read_csv(str(PROJECT_ROOT / paths_config["partial_metadata"]))

MERGED_CSV = PROJECT_ROOT / paths_config["merged_csv"]



GII_3 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_vp1"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["vp1_gii_3"],
    "genotype",
    "GII.3"
)

GII_4 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_vp1"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["vp1_gii_4"],
    "genotype",
    "GII.4"
)

GII_6 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_vp1"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["vp1_gii_6"],
    "genotype",
    "GII.6"
)

# ######################

## runs cd hit on common genotypes/p-types to cluster
CLUSTERED_GII_3 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["vp1_gii_3"],
    PROJECT_ROOT / paths_config["vp1_clustered_gii_3"],
    "0.95"
)

CLUSTERED_GII_4 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["vp1_gii_4"],
    PROJECT_ROOT / paths_config["vp1_clustered_gii_4"],
    "0.95"
)

CLUSTERED_GII_6 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["vp1_gii_6"],
    PROJECT_ROOT / paths_config["vp1_clustered_gii_6"],
)

## set containing all accessions from complete sequences
raw_set_vp1 = cd.get_accesions(PROJECT_ROOT / paths_config["final_vp1"])



#########################################
VP1_BEFORE_CDHIT_TSV = PROJECT_ROOT / paths_config["vp1_before_cdhit_tsv"]

with open(VP1_BEFORE_CDHIT_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['accession'])
        for accession in raw_set_vp1:
            writer.writerow([accession])


#########################################




# ## set for all accesions for sequences that were sent to cd hit:
common_types_accessions_vp1 = cd.get_accesions(
    PROJECT_ROOT / paths_config["vp1_gii_3"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["vp1_gii_4"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["vp1_gii_6"]
)
 

# ## retrive accession for sequences not sent to cd hit
rare_types_accessions_vp1 = raw_set_vp1.difference(common_types_accessions_vp1)


## retrieve set for accession of all clustered sequences:
clustered_accessions_vp1 = cd.get_accesions(
    PROJECT_ROOT / paths_config["vp1_clustered_gii_3"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["vp1_clustered_gii_4"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["vp1_clustered_gii_6"]
)

## set for accesions that are used in the final global alignment
final_vp1_accessions = rare_types_accessions_vp1 | clustered_accessions_vp1
print('------------------------')
print(f'antal sekvenser från complete (som använts i global msa) = {len(final_vp1_accessions)}')

FINAL_VP1_ACCESSION_TSV = PROJECT_ROOT / paths_config["final_vp1_accession_tsv"]

with open(FINAL_VP1_ACCESSION_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['accession'])
        for accession in final_vp1_accessions:
            writer.writerow([accession])


VP1 = cd.get_sequences(
    PROJECT_ROOT / paths_config["final_vp1"], 
    PROJECT_ROOT / paths_config["vp1"],
    final_vp1_accessions
)
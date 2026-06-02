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


# CLUSTERED_RDRP = cd.cluster_sequences(
#     PROJECT_ROOT / paths_config["final_rdrp"], 
#     PROJECT_ROOT / paths_config["clustered_rdrp"])



df_complete = pd.read_csv(str(PROJECT_ROOT / paths_config["input_metadata"]))
df_partial = pd.read_csv(str(PROJECT_ROOT / paths_config["partial_metadata"]))
MERGED_CSV = PROJECT_ROOT / paths_config["merged_csv"]
merged = pd.concat([df_complete, df_partial], ignore_index=True)

merged = merged.drop_duplicates(subset = "accession")
merged.to_csv(MERGED_CSV, index = False)

GII_P7 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_rdrp"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["rdrp_gii_p7"],
    "p_type",
    "GII.P7"
)

GII_P16 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_rdrp"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["rdrp_gii_p16"],
    "p_type",
    "GII.P16"
)

GII_P21 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_rdrp"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["rdrp_gii_p21"],
    "p_type",
    "GII.P21"
)

GII_P31 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["final_rdrp"],
    MERGED_CSV,
    PROJECT_ROOT / paths_config["rdrp_gii_p31"],
    "p_type",
    "GII.P31"
)

######################

## runs cd hit on common genotypes/p-types to cluster
CLUSTERED_GII_P7 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["rdrp_gii_p7"],
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p7"],
    "0.9725"
)

CLUSTERED_GII_P16 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["rdrp_gii_p16"],
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p16"],
    "0.965"
)
CLUSTERED_GII_P21 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["rdrp_gii_p21"],
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p21"],
    "0.98"
)

CLUSTERED_GII_P31 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["rdrp_gii_p31"],
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p31"],
    "0.975"
)






## set containing all accessions from complete sequences
raw_set = cd.get_accesions(PROJECT_ROOT / paths_config["final_rdrp"])


## set for all accesions for sequences that were sent to cd hit:
common_types_accessions = cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_gii_p7"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_gii_p16"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_gii_p21"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_gii_p31"]
) 
 

## retrive accession for sequences not sent to cd hit
rare_types_accessions = raw_set.difference(common_types_accessions)


## retrieve set for accession of all clustered sequences:
clustered_accessions = cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p7"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p16"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p21"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["rdrp_clustered_gii_p31"]
) 

## set for accesions that are used in the final global alignment
final_rdrp_accessions = rare_types_accessions | clustered_accessions
print('------------------------')
print(f'antal sekvenser från complete (som använts i global msa) = {len(final_rdrp_accessions)}')

FINAL_RDRP_ACCESSION_TSV = PROJECT_ROOT / paths_config["final_rdrp_accession_tsv"]

with open(FINAL_RDRP_ACCESSION_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        writer.writerow(['accession'])
        for accession in final_rdrp_accessions:
            writer.writerow([accession])


POLYMERASE = cd.get_sequences(
    PROJECT_ROOT / paths_config["final_rdrp"], 
    PROJECT_ROOT / paths_config["polymerase"],
    final_rdrp_accessions
)
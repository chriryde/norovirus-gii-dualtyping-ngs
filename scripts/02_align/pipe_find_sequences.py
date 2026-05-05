import alignment_tools as at
import cd_hit as cd
import partial as pa
from pathlib import Path
import yaml
import csv

ALIGN_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (ALIGN_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "alignment.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

paths_config = config["paths"]


# ## get rdrp region
RDRP_FASTA = at.get_region(
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT /paths_config["input_fasta"],
    PROJECT_ROOT /paths_config["rdrp_fasta"],
    "rdrp_start",
    "rdrp_end"
)

 ## get vp1 region
VP1_FASTA = at.get_region(
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["input_fasta"],
    PROJECT_ROOT / paths_config["vp1_fasta"],
    "vp1_start",
    "vp1_end"
)

# ## get rdrp+vp1 region
RDRP_VP1_FASTA = at.get_region(
    PROJECT_ROOT /paths_config["input_metadata"],
    PROJECT_ROOT /paths_config["input_fasta"],
    PROJECT_ROOT /paths_config["rdrp_vp1_fasta"],
    "rdrp_start",
    "vp1_end"
)


# partial_accessions = pa.partial_accessions(
#     PROJECT_ROOT / paths_config["input_metadata"],
#     PROJECT_ROOT / paths_config["input_metadata"]
# )

# PARTIAL_FASTA = cd.get_sequences(
#     PROJECT_ROOT / paths_config["input_partial"],
#     PROJECT_ROOT / paths_config["partial_fasta"],
#     partial_accessions
# )

# FINAL_RDRP = pa.m(
    
# )

###################################
# ## runs get genotype
GII_4 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_4"],
    "genotype",
    "GII.4"
)

GII_3 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_3"],
    "genotype",
    "GII.3"
)

GII_2 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_2"],
    "genotype",
    "GII.2"
)

GII_17 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_17"],
    "genotype",
    "GII.17"
)

# # ## runs get p_type

GII_P4 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p4"],
    "p_type",
    "GII.P4"
)

GII_P16 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p16"],
    "p_type",
    "GII.P16"
)


GII_P17 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p17"],
    "p_type",
    "GII.P17"
)

GII_P21 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p21"],
    "p_type",
    "GII.P21"
)

GII_P31 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p31"],
    "p_type",
    "GII.P31"
)



## run cd hit to cluster
CLUSTERED_GII_2 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_2"],
    PROJECT_ROOT / paths_config["clustered_gii_2"]
)

CLUSTERED_GII_3 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_3"],
    PROJECT_ROOT / paths_config["clustered_gii_3"]
)

CLUSTERED_GII_4 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_4"],
    PROJECT_ROOT / paths_config["clustered_gii_4"],
    "0.95"
)

CLUSTERED_GII_17 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_17"],
    PROJECT_ROOT / paths_config["clustered_gii_17"]
)


CLUSTERED_GII_P4 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_p16"],
    PROJECT_ROOT / paths_config["clustered_gii_p4"]
)

CLUSTERED_GII_P16 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_p16"],
    PROJECT_ROOT / paths_config["clustered_gii_p16"]
)



CLUSTERED_GII_P17 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_p17"],
    PROJECT_ROOT / paths_config["clustered_gii_p17"]
)


CLUSTERED_GII_P21 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_p21"],
    PROJECT_ROOT / paths_config["clustered_gii_p21"]
)


CLUSTERED_GII_P16 = cd.cluster_sequences(
    PROJECT_ROOT / paths_config["gii_p31"],
    PROJECT_ROOT / paths_config["clustered_gii_p31"]
)


## set containing all accessions
raw_set = cd.get_accesions(PROJECT_ROOT / paths_config["input_fasta"])


## set for all accesions for sequences that were sent to cd hit:
common_types_accessions = cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_2"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_3"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_4"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_17"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_p4"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_p16"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_p17"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_p21"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["gii_p31"]
) 
 

## retrive accession for sequences not sent to cd hit
rare_types_accessions = raw_set.difference(common_types_accessions)


## retrieve set for accession of all clustered sequences:
clustered_accessions = cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_2"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_3"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_4"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_17"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_p4"]
)  | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_p16"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_p17"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_p21"]
) | cd.get_accesions(
    PROJECT_ROOT / paths_config["clustered_gii_p31"]
) 


final_accesions = rare_types_accessions | clustered_accessions

ACCESSION_TSV = PROJECT_ROOT / paths_config["accession_tsv"]

with open(ACCESSION_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for accession in final_accesions:
            writer.writerow([accession])

# print(f'raw = {len(raw_set)}')
# print(f'common = {len(common_types_accessions)}')
# print(f'rare = {len(rare_types_accessions)}')
# print(f'clustered = {len(clustered_accessions)}')
# print(f'final = {len(final_accesions)}')


## create new fasta file with all rare types + clustered sequences
FINAL_GLOBAL = cd.get_sequences(
    PROJECT_ROOT / paths_config["input_fasta"], 
    PROJECT_ROOT / paths_config["final_global"],
    final_accesions
)

FINAL_RDRP = cd.get_sequences(
    PROJECT_ROOT / paths_config["rdrp_fasta"], 
    PROJECT_ROOT / paths_config["final_rdrp"],
    final_accesions
)

FINAL_VP1 = cd.get_sequences(
    PROJECT_ROOT / paths_config["vp1_fasta"], 
    PROJECT_ROOT / paths_config["final_vp1"],
    final_accesions
)


FINAL_RDRP_VP1 = cd.get_sequences(
    PROJECT_ROOT / paths_config["rdrp_vp1_fasta"], 
    PROJECT_ROOT / paths_config["final_rdrp_vp1"],
    final_accesions
)
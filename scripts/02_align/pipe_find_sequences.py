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
RDRP_COMPLETE_FASTA = at.get_region(
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT /paths_config["input_fasta"],
    PROJECT_ROOT /paths_config["rdrp_complete_fasta"],
    "rdrp_start",
    "rdrp_end"
)

 ## get vp1 region
VP1_COMPLETE_FASTA = at.get_region(
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["input_fasta"],
    PROJECT_ROOT / paths_config["vp1_complete_fasta"],
    "vp1_start",
    "vp1_end"
)

# ## get rdrp+vp1 region
# RDRP_VP1_FASTA = at.get_region(
#     PROJECT_ROOT /paths_config["input_metadata"],
#     PROJECT_ROOT /paths_config["input_fasta"],
#     PROJECT_ROOT /paths_config["rdrp_vp1_fasta"],
#     "rdrp_start",
#     "vp1_end"
# )

#### 
partial_rdrp, partial_vp1 = pa.partial_accessions(
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["partial_metadata"]
)

## alla sekvenser från partial som vi vill komplettera med
# PARTIAL_FASTA = cd.get_sequences(
#     PROJECT_ROOT / paths_config["input_partial"],
#     PROJECT_ROOT / paths_config["partial_fasta"],
#     partial_rdrp | partial_vp1
# )

## rdrp regioner från utvalda partial
PARTIAL_RDRP_FASTA = cd.get_sequences(
    PROJECT_ROOT / paths_config["input_partial"],
    PROJECT_ROOT / paths_config["partial_rdrp_fasta"],
    partial_rdrp
)

## vp1 regioner från utvalda partials
PARTIAL_VP1_FASTA = cd.get_sequences(
    PROJECT_ROOT / paths_config["input_partial"],
    PROJECT_ROOT / paths_config["partial_vp1_fasta"],
    partial_vp1
)

###################################

## plockar ut alla complete sekvenser inom de mest förekommande genotyperna
# ## runs get genotype
GII_4 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_4"],
    "genotype",
    "GII.4"
)

GII_3 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_3"],
    "genotype",
    "GII.3"
)

GII_2 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_2"],
    "genotype",
    "GII.2"
)

GII_17 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["vp1_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_17"],
    "genotype",
    "GII.17"
)

# # ## runs get p_type
## plockar ut alla complete sekvenser från de mest förekommande p-typerna
GII_P4 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p4"],
    "p_type",
    "GII.P4"
)

GII_P16 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p16"],
    "p_type",
    "GII.P16"
)


GII_P17 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p17"],
    "p_type",
    "GII.P17"
)

GII_P21 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p21"],
    "p_type",
    "GII.P21"
)

GII_P31 = at.get_genotype_or_ptype(
    PROJECT_ROOT / paths_config["rdrp_complete_fasta"],
    PROJECT_ROOT / paths_config["input_metadata"],
    PROJECT_ROOT / paths_config["gii_p31"],
    "p_type",
    "GII.P31"
)

######################

## runs cd hit on common genotypes/p-types to cluster
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


## set containing all accessions from complete sequences
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

## set for accesions that are used in the final global alignment
final_accesions = rare_types_accessions | clustered_accessions
print('------------------------')
print(f'antal sekvenser från complete (som använts i global msa) = {len(final_accesions)}')

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
## to be used in the global msa
FINAL_GLOBAL = cd.get_sequences(
    PROJECT_ROOT / paths_config["input_fasta"], 
    PROJECT_ROOT / paths_config["final_global"],
    final_accesions
)


#### clustered rdrp and vp1

## rdrp regionr från sekvenser i global msa
FILTERED_RDRP = cd.get_sequences(
    PROJECT_ROOT / paths_config["rdrp_complete_fasta"], 
    PROJECT_ROOT / paths_config["filtered_rdrp"],
    final_accesions
)

## vp1 regioner från sekvenser i global msa
FILTERED_VP1 = cd.get_sequences(
    PROJECT_ROOT / paths_config["vp1_complete_fasta"], 
    PROJECT_ROOT / paths_config["filtered_vp1"],
    final_accesions
)


# FINAL_RDRP_VP1 = cd.get_sequences(
#     PROJECT_ROOT / paths_config["rdrp_vp1_fasta"], 
#     PROJECT_ROOT / paths_config["final_rdrp_vp1"],
#     final_accesions
# )

############ hämtar rätt regioner från partial
### rdrp regioner från partial sekvenser
RDRP_PARTIAL = at.get_region(
    PROJECT_ROOT / paths_config["partial_metadata"],
    PROJECT_ROOT /paths_config["partial_rdrp_fasta"],
    PROJECT_ROOT /paths_config["rdrp_partial"],
    "rdrp_start",
    "rdrp_end"
)

### vp1 regioner från partial sekvenser
VP1_PARTIAL = at.get_region(
    PROJECT_ROOT / paths_config["partial_metadata"],
    PROJECT_ROOT /paths_config["partial_vp1_fasta"],
    PROJECT_ROOT /paths_config["vp1_partial"],
    "vp1_start",
    "vp1_end"
)


# ## merge rdrp regioner från complete & partial sekvenser
# ## fil som ska användas i msa
# FINAL_RDRP = pa.merge_partial_regions(
#     RDRP_PARTIAL,
#     FILTERED_RDRP,
#     PROJECT_ROOT / paths_config["final_rdrp"]
# )

# ## merge vp1 regioner från complete & partial sekvenser
# ## fil som ska användas i msa
# FINAL_VP1 = pa.merge_partial_regions(
#     PARTIAL_FASTA,
#     FILTERED_VP1,
#     PROJECT_ROOT / paths_config["final_vp1"]
# )

## merge rdrp regioner från complete & partial sekvenser
## fil som ska användas i msa
FINAL_RDRP = pa.merge_partial_regions(
    PARTIAL_RDRP_FASTA,
    FILTERED_RDRP,
    PROJECT_ROOT / paths_config["final_rdrp"]
)

## merge vp1 regioner från complete & partial sekvenser
## fil som ska användas i msa
FINAL_VP1 = pa.merge_partial_regions(
    PARTIAL_VP1_FASTA,
    FILTERED_VP1,
    PROJECT_ROOT / paths_config["final_vp1"]
)
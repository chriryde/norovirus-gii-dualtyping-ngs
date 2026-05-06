import curation_tools as ct
import pandas as pd
import yaml

from pathlib import Path
from Bio import SeqIO

def project_path(pathlike: str | Path) -> Path:
    path = Path(pathlike)
    return path if path.is_absolute() else PROJECT_ROOT / path

CURATE_DIR: Path = Path(__file__).resolve().parent
PROJECT_ROOT: Path = (CURATE_DIR.parent).parent
CONFIG_PATH: Path = PROJECT_ROOT / "config" / "curation.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

METADATA_DIR = project_path(config["paths"]["metadata_dir"])

INPUT_METADATA = project_path(config["paths"]["input_metadata"])
RAW_FASTA = project_path(config["paths"]["input_fasta"])
INPUT_CDS_FASTA = project_path(config["paths"]["input_cds_fasta"])
TEMP_TYPING_DIR = project_path(config["paths"]["temp_typing_dir"])

COMPLETE_FASTA_PATH = project_path(config["paths"]["complete_fasta"])
COMPLETE_FILTERED_FASTA_PATH = project_path(config["paths"]["complete_filtered_fasta"])
COMPLETE_CDS_FASTA_PATH = project_path(config["paths"]["complete_cds_fasta"])

PARTIAL_FASTA_PATH = project_path(config["paths"]["partial_fasta"])
PARTIAL_FILTERED_FASTA_PATH = project_path(config["paths"]["partial_filtered_fasta"])
PARTIAL_CDS_FASTA_PATH = project_path(config["paths"]["partial_cds_fasta"])

HUCAT_GENBANK = project_path(config["paths"]["hucat_genbank"])
HUCAT_FASTA = project_path(config["paths"]["hucat_fasta"])

HUCAT_COMPLETE_FASTA_PATH = project_path(config["paths"]["hucat_complete_fasta"])
HUCAT_COMPLETE_FILTERED_FASTA_PATH = project_path(config["paths"]["hucat_complete_filtered_fasta"])
HUCAT_COMPLETE_FILTERED_GENBANK_PATH = project_path(config["paths"]["hucat_complete_filtered_genbank"])

HUCAT_PARTIAL_FASTA_PATH = project_path(config["paths"]["hucat_partial_fasta"])
HUCAT_PARTIAL_FILTERED_FASTA_PATH = project_path(config["paths"]["hucat_partial_filtered_fasta"])
HUCAT_PARTIAL_FILTERED_GENBANK_PATH = project_path(config["paths"]["hucat_partial_filtered_genbank"])

FINAL_COMPLETE_FASTA = project_path(config["paths"]["final_complete_fasta"])
FINAL_PARTIAL_FASTA = project_path(config["paths"]["final_partial_fasta"])

if not TEMP_TYPING_DIR.is_dir():
    print("shouldn't run this")

    complete_sequences, partial_sequences = ct.get_sequence_completeness(
        INPUT_METADATA
    )
    print(len(complete_sequences))

    SeqIO.convert(HUCAT_GENBANK, "genbank", HUCAT_FASTA, "fasta")

    genbank_complete_sequences, genbank_partial_sequences = ct.get_genbank_sequence_completeness(
        HUCAT_GENBANK
    )

    partial_sequences = partial_sequences - genbank_complete_sequences
    print(len(partial_sequences))

    hucat_complete_sequences = genbank_complete_sequences - complete_sequences
    hucat_partial_sequences = genbank_partial_sequences - partial_sequences

# #|===| START OF COMPLETE BLOCK |===============================================|
    COMPLETE_FASTA: Path = ct.extract_and_save_to_fasta(
        RAW_FASTA, 
        complete_sequences,
        COMPLETE_FASTA_PATH
    )

    complete_gii_sequences: set[str] = ct.get_GII_sequences(
        INPUT_METADATA,
        complete_sequences
    )
    assert(complete_gii_sequences)

    complete_length_filtered_sequences: set[str] = ct.get_length_filtered_sequences(
        COMPLETE_FASTA,
        config["criteria"]["length_threshold"]
    )
    assert(complete_length_filtered_sequences)

    complete_unique_sequences: set[str] = ct.get_unique_sequences(
        COMPLETE_FASTA
    )
    assert(complete_unique_sequences)

    complete_annotated_sequences: set[str] = ct.get_annotated_sequences(
        INPUT_METADATA
    )
    assert(complete_annotated_sequences)

    complete_non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        COMPLETE_FASTA,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(complete_non_ambiguous_sequences)

    complete_filtered_sequences: set[str] = complete_sequences \
        & complete_length_filtered_sequences & complete_unique_sequences \
        & complete_annotated_sequences & complete_non_ambiguous_sequences \
        & complete_gii_sequences
    assert complete_filtered_sequences

    COMPLETE_FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        COMPLETE_FASTA, 
        complete_filtered_sequences,
        COMPLETE_FILTERED_FASTA_PATH
    )

    excluded_accessions = complete_sequences - complete_filtered_sequences
    ct.export_excluded_sequences(
        project_path(config["paths"]["excluded_csv"]),
        excluded_accessions=excluded_accessions,
        gii_sequences=complete_gii_sequences,
        complete_sequences=complete_sequences,
        length_filtered_sequences=complete_length_filtered_sequences,
        unique_sequences=complete_unique_sequences,
        annotated_sequences=complete_annotated_sequences,
        non_ambiguous_sequences=complete_non_ambiguous_sequences
    )

    # ct.typing_tool_intialise(
    #     COMPLETE_FILTERED_FASTA,
    #     "complete",
    #     batch_size=500
    # )

    COMPLETE_CDS_FASTA: Path = ct.extract_and_save_to_fasta(
        INPUT_CDS_FASTA, 
        complete_filtered_sequences,
        COMPLETE_CDS_FASTA_PATH
    )
#|===| END OF COMPLETE BLOCK |================================================|
#|                                                                            |
#|===| START OF PARTIAL BLOCK |===============================================|
    PARTIAL_FASTA: Path = ct.extract_and_save_to_fasta(
        RAW_FASTA, 
        partial_sequences,
        PARTIAL_FASTA_PATH
    )

    PARTIAL_CDS_FASTA: Path = ct.extract_and_save_to_fasta(
        INPUT_CDS_FASTA, 
        partial_sequences,
        PARTIAL_CDS_FASTA_PATH
    )

    PARTIAL_REGION_INFORMATION: Path | None = ct.get_genomic_info(
        project_path(config["paths"]["partial_region_information"]),
        typing_information=False,
        region_information=True,
        CDS_FASTA=PARTIAL_CDS_FASTA
    )
    if PARTIAL_REGION_INFORMATION is None:
        raise RuntimeError("Failed to create partial region information")

    partials_rdrp: set[str] = ct.get_rdrp_sequences(
        PARTIAL_REGION_INFORMATION,
    )

    partials_vp1: set[str] = ct.get_vp1_sequences(
        PARTIAL_REGION_INFORMATION,
    )

    partials_junction: set[str] = ct.get_junction_sequences(
        PARTIAL_REGION_INFORMATION,
    )

    partial_gii_sequences: set[str] = ct.get_GII_sequences(
        project_path(config["paths"]["input_metadata"]),
        partial_sequences
    )
    assert(partial_gii_sequences)

    partial_unique_sequences: set[str] = ct.get_unique_sequences(
        PARTIAL_FASTA_PATH
    )
    assert(partial_unique_sequences)

    partial_non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        PARTIAL_FASTA_PATH,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(partial_non_ambiguous_sequences)

    partial_rdrp_filtered_sequences: set[str] = partial_sequences \
        & partials_rdrp & partial_gii_sequences & partial_unique_sequences \
        & partial_non_ambiguous_sequences
    assert(partial_rdrp_filtered_sequences)

    partial_vp1_filtered_sequences: set[str] = partial_sequences \
        & partials_vp1 & partial_gii_sequences & partial_unique_sequences \
        & partial_non_ambiguous_sequences
    assert(partial_vp1_filtered_sequences)

    partial_junction_filtered_sequences: set[str] = partial_sequences \
        & partials_junction & partial_gii_sequences & partial_unique_sequences \
        & partial_non_ambiguous_sequences
    assert(partial_junction_filtered_sequences)

    partial_all_filtered_sequences = partial_rdrp_filtered_sequences \
        | partial_vp1_filtered_sequences | partial_junction_filtered_sequences

    PARTIAL_FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        PARTIAL_FASTA_PATH, 
        partial_all_filtered_sequences,
        PARTIAL_FILTERED_FASTA_PATH
    )

    # ct.typing_tool_intialise(
    #     PARTIAL_FILTERED_FASTA,
    #     "partial",
    #     batch_size=500
    # )
#|===| END OF PARTIAL BLOCK |=================================================|

#|===| START OF COMPLETE GENBANK BLOCK |======================================|
    HUCAT_COMPLETE_FASTA: Path = ct.extract_and_save_to_fasta(
        HUCAT_FASTA, 
        hucat_complete_sequences,
        HUCAT_COMPLETE_FASTA_PATH
    )

    hucat_complete_gii_sequences: set[str] = ct.get_GII_sequences(
        HUCAT_GENBANK,
        hucat_complete_sequences,
        genbank_mode=True
    )
    assert(hucat_complete_gii_sequences)

    hucat_complete_length_filtered_sequences: set[str] = ct.get_length_filtered_sequences(
        HUCAT_COMPLETE_FASTA,
        config["criteria"]["length_threshold"]
    )
    assert(hucat_complete_length_filtered_sequences)

    hucat_complete_unique_sequences: set[str] = ct.get_unique_sequences(
        HUCAT_COMPLETE_FASTA
    )
    assert(hucat_complete_unique_sequences)

    hucat_complete_non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        HUCAT_COMPLETE_FASTA,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(hucat_complete_non_ambiguous_sequences)
    
    hucat_complete_filtered_sequences: set[str] = hucat_complete_sequences \
        & hucat_complete_length_filtered_sequences & hucat_complete_unique_sequences \
        & hucat_complete_non_ambiguous_sequences & hucat_complete_gii_sequences
    assert hucat_complete_filtered_sequences

    HUCAT_COMPLETE_FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        HUCAT_COMPLETE_FASTA, 
        hucat_complete_filtered_sequences,
        HUCAT_COMPLETE_FILTERED_FASTA_PATH
    )

    # ct.typing_tool_intialise(
    #     HUCAT_COMPLETE_FILTERED_FASTA,
    #     "hucat_complete",
    #     batch_size=500
    # )

    HUCAT_COMPLETE_FILTERED_GENBANK: Path = ct.extract_and_save_to_genbank(
        HUCAT_GENBANK, 
        hucat_complete_filtered_sequences,
        HUCAT_COMPLETE_FILTERED_GENBANK_PATH
    )
#|===| END OF COMPLETE GENBANK BLOCK |========================================|
#|                                                                            |
#|===| START OF PARTIAL GENBANK BLOCK |=======================================|
    HUCAT_PARTIAL_FASTA: Path = ct.extract_and_save_to_fasta(
        HUCAT_FASTA, 
        hucat_partial_sequences,
        HUCAT_PARTIAL_FASTA_PATH
    )

    hucat_partial_gii_sequences: set[str] = ct.get_GII_sequences(
        HUCAT_GENBANK,
        hucat_partial_sequences,
        genbank_mode=True
    )
    assert(hucat_partial_gii_sequences)

    hucat_partial_unique_sequences: set[str] = ct.get_unique_sequences(
        HUCAT_PARTIAL_FASTA
    )
    assert(hucat_partial_unique_sequences)

    hucat_partial_non_ambiguous_sequences: set[str] = ct.get_ambiguous_filtered_sequences(
        HUCAT_PARTIAL_FASTA,
        config["criteria"]["ambiguous_threshold"]
    )
    assert(hucat_partial_non_ambiguous_sequences)

    HUCAT_PARTIAL_REGION_INFORMATION: Path | None = ct.get_genomic_info(
        project_path(config["paths"]["hucat_partial_region_information"]),
        typing_information=False,
        region_information=True,
        genbank_mode=True,
        GENBANK_FILE=HUCAT_GENBANK
    )
    if HUCAT_PARTIAL_REGION_INFORMATION is None:
        raise RuntimeError("Failed to create partial region information")
    
    hucat_partials_rdrp: set[str] = ct.get_rdrp_sequences(
        HUCAT_PARTIAL_REGION_INFORMATION,
    )
    assert(hucat_partials_rdrp)

    hucat_partials_vp1: set[str] = ct.get_vp1_sequences(
        HUCAT_PARTIAL_REGION_INFORMATION,
    )
    assert(hucat_partials_vp1)

    hucat_partials_junction: set[str] = ct.get_junction_sequences(
        HUCAT_PARTIAL_REGION_INFORMATION,
    )
    assert(hucat_partials_junction)

    hucat_partial_rdrp_filtered_sequences: set[str] = hucat_partial_sequences \
        & hucat_partials_rdrp & hucat_partial_unique_sequences \
        & hucat_partial_non_ambiguous_sequences & hucat_partial_gii_sequences
    assert(hucat_partial_rdrp_filtered_sequences)

    hucat_partial_vp1_filtered_sequences: set[str] = hucat_partial_sequences \
        & hucat_partials_vp1 & hucat_partial_unique_sequences \
        & hucat_partial_non_ambiguous_sequences & hucat_partial_gii_sequences
    assert(hucat_partial_vp1_filtered_sequences)

    hucat_partial_junction_filtered_sequences: set[str] = hucat_partial_sequences \
        & hucat_partials_junction & hucat_partial_unique_sequences \
        & hucat_partial_non_ambiguous_sequences & hucat_partial_gii_sequences
    assert(hucat_partial_junction_filtered_sequences)

    hucat_partial_all_filtered_sequences = hucat_partial_rdrp_filtered_sequences \
        | hucat_partial_vp1_filtered_sequences | hucat_partial_junction_filtered_sequences
    
    HUCAT_PARTIAL_FILTERED_FASTA: Path = ct.extract_and_save_to_fasta(
        HUCAT_PARTIAL_FASTA, 
        hucat_partial_all_filtered_sequences,
        HUCAT_PARTIAL_FILTERED_FASTA_PATH
    )

    # ct.typing_tool_intialise(
    #     HUCAT_PARTIAL_FILTERED_FASTA,
    #     "hucat_partial",
    #     batch_size=500
    # )

    HUCAT_PARTIAL_FILTERED_GENBANK: Path = ct.extract_and_save_to_genbank(
        HUCAT_GENBANK, 
        hucat_partial_all_filtered_sequences,
        HUCAT_PARTIAL_FILTERED_GENBANK_PATH
    )
#|===| END OF PARTIAL GENBANK BLOCK |=========================================|
#|                                                                            |
#|===| START OF COMPLETE BLOCK |==============================================|
COMPLETE_TYPING_CSV: Path | None = ct.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "complete"
)
if COMPLETE_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

COMPLETE_TYPING_CSV = project_path(COMPLETE_TYPING_CSV)

typing_df = pd.read_csv(COMPLETE_TYPING_CSV)

COMPLETE_GENOMIC_SPECIFICATION: Path | None = ct.get_genomic_info(
    project_path(config["paths"]["complete_full_specification"]),
    typing_information=True,
    TYPING_FILE=COMPLETE_TYPING_CSV,
    region_information=True,
    CDS_FASTA=COMPLETE_CDS_FASTA_PATH,
)
# #|===| END OF COMPLETE BLOCK |================================================|
# #|                                                                            |
# #|===| START OF PARTIAL BLOCK |===============================================|
PARTIAL_CDS_FASTA = project_path(config["paths"]["partial_cds_fasta"])

PARTIAL_TYPING_CSV: Path | None = ct.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "partial"
)
if PARTIAL_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

PARTIAL_TYPING_CSV = project_path(PARTIAL_TYPING_CSV)

PARTIAL_SPECIFICATION: Path | None = ct.get_genomic_info(
    project_path(config["paths"]["partial_full_specification"]),
    typing_information=True,
    TYPING_FILE=PARTIAL_TYPING_CSV,
    region_information=True,
    CDS_FASTA=PARTIAL_CDS_FASTA,
)
#|===| END OF PARTIAL BLOCK |=================================================|
#|                                                                            |
#|===| START OF COMPLETE GENBANK BLOCK |======================================|
HUCAT_COMPLETE_TYPING_CSV: Path | None = ct.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "hucat_complete"
)
if HUCAT_COMPLETE_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

HUCAT_COMPLETE_TYPING_CSV = project_path(HUCAT_COMPLETE_TYPING_CSV)

HUCAT_COMPLETE_SPECIFICATION: Path | None = ct.get_genomic_info(
    project_path(config["paths"]["hucat_complete_full_specification"]),
    typing_information=True,
    TYPING_FILE=HUCAT_COMPLETE_TYPING_CSV,
    region_information=True,
    genbank_mode=True,
    GENBANK_FILE=HUCAT_COMPLETE_FILTERED_GENBANK_PATH
)
#|===| END OF COMPLETE GENBANK BLOCK |========================================|
#|                                                                            |
#|===| START OF PARTIAL GENBANK BLOCK |=======================================|
HUCAT_PARTIAL_TYPING_CSV: Path | None = ct.typing_tool_get_results(
    TEMP_TYPING_DIR,
    "hucat_partial"
)
if HUCAT_PARTIAL_TYPING_CSV is None:
    print(f"Typing tool results are not finished yet. Halting pipeline.")
    raise SystemExit

HUCAT_PARTIAL_TYPING_CSV = project_path(HUCAT_PARTIAL_TYPING_CSV)

HUCAT_PARTIAL_SPECIFICATION: Path | None = ct.get_genomic_info(
    project_path(config["paths"]["hucat_partial_full_specification"]),
    typing_information=True,
    TYPING_FILE=HUCAT_PARTIAL_TYPING_CSV,
    region_information=True,
    genbank_mode=True,
    GENBANK_FILE=HUCAT_PARTIAL_FILTERED_GENBANK_PATH
)
#|===| END OF PARTIAL GENBANK BLOCK |=========================================|

FULL_COMPLETE_SPECIFICATION = project_path(config["paths"]["final_complete_full_specification"])

complete_csv_files = [
    project_path(config["paths"]["complete_full_specification"]),
    project_path(config["paths"]["hucat_complete_full_specification"])
]

complete_df = pd.concat(
    [pd.read_csv(file) for file in complete_csv_files],
    ignore_index=True
)

complete_df = ct.prepare_final_complete_specification(complete_df)

print(
    complete_df[complete_df["accession"] == "AB541286.1"]
    .to_dict("records")
)

complete_df.to_csv(FULL_COMPLETE_SPECIFICATION, index=False)

final_complete_sequences = set(
    complete_df["accession"].dropna().astype(str).str.strip()
)

FULL_PARTIAL_SPECIFICATION = project_path(config["paths"]["final_partial_full_specification"])

partial_csv_files = [
    project_path(config["paths"]["partial_full_specification"]),
    project_path(config["paths"]["hucat_partial_full_specification"]),
]

partial_df = pd.concat(
    [pd.read_csv(file) for file in partial_csv_files],
    ignore_index=True
)

partial_df.to_csv(FULL_PARTIAL_SPECIFICATION, index=False)
FULL_PARTIAL_SPECIFICATION = ct.filter_csv_on_ptype_genotype(FULL_PARTIAL_SPECIFICATION, FULL_PARTIAL_SPECIFICATION)
partial_df = pd.read_csv(FULL_PARTIAL_SPECIFICATION)

final_partial_sequences = set(partial_df["accession"].dropna().astype(str).str.strip())



MERGED_COMPLETE_FASTA = project_path(config["paths"]["merged_complete_fasta"])
MERGED_PARTIAL_FASTA = project_path(config["paths"]["merged_partial_fasta"])

ct.merge_fasta(
    [COMPLETE_FILTERED_FASTA_PATH, HUCAT_COMPLETE_FILTERED_FASTA_PATH],
    MERGED_COMPLETE_FASTA
)

ct.merge_fasta(
    [PARTIAL_FILTERED_FASTA_PATH, HUCAT_PARTIAL_FILTERED_FASTA_PATH],
    MERGED_PARTIAL_FASTA
)

ct.extract_and_save_to_fasta(
    MERGED_COMPLETE_FASTA,
    final_complete_sequences,
    FINAL_COMPLETE_FASTA
)

ct.extract_and_save_to_fasta(
    MERGED_PARTIAL_FASTA,
    final_partial_sequences,
    FINAL_PARTIAL_FASTA
)

final_partial_sequences = ct.get_accessions(FINAL_PARTIAL_FASTA)

FULL_PARTIAL_SPECIFICATION = ct.filter_csv_accessions(
    FULL_PARTIAL_SPECIFICATION, 
    FULL_PARTIAL_SPECIFICATION,
    final_partial_sequences
)
partial_df = pd.read_csv(FULL_PARTIAL_SPECIFICATION)

complete_ass = ct.assert_accessions(FINAL_COMPLETE_FASTA, FULL_COMPLETE_SPECIFICATION)
print(len(complete_ass))

partial_ass = ct.assert_accessions(FINAL_PARTIAL_FASTA, FULL_PARTIAL_SPECIFICATION)
print(len(partial_ass))

# ct.clean_up_temporary_files(TYPING_CSV)
import pandas as pd
import subprocess
import datetime
import shutil
import json
import math
import time
import csv

from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError
from Bio.SeqFeature import BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import Any
from Bio import SeqIO

vp1_names: set[str] = {
    "VP1",
    "vp1",
    "capsid VP1",
    "viral protein 1",
    "capsid protein VP1",
    "major capsid protein",
    "major viral capsid protein",
    "major capsid protein VP1"
}

rdrp_names: set[str] = {
    "RdRp",
    "rdrp",
    "RNA-dependent RNA polymerase"
}

polyprotein_names: set[str] = {
    "polyprotein",
    "nonstructural polyprotein"
}

def make_row_full_information():
    return {
        "BLAST_score"     : None,
        "begin"           : None,
        "end"             : None,
        "length"          : None,
        "p_type"          : None,
        "p_subtype"       : None,
        "rdrp_start"      : None,
        "rdrp_end"        : None,
        "rdrp_source"     : None,
        "genotype"        : None,
        "genotype_subtype": None,
        "vp1_start"       : None,
        "vp1_end"         : None,
        "has_overlap"     : None,
        "overlap"         : None
    }

def make_row_region_information():
    return {
        "rdrp_start"      : None,
        "rdrp_end"        : None,
        "rdrp_source"     : None,
        "vp1_start"       : None,
        "vp1_end"         : None,
        "has_overlap"     : None,
        "overlap"         : None
    }

def make_row_typing_information():
    return {
        "BLAST_score"     : None,
        "begin"           : None,
        "end"             : None,
        "length"          : None,
        "p_type"          : None,
        "p_subtype"       : None,
        "genotype"        : None,
        "genotype_subtype": None
    }


def get_accessions(INPUT_FASTA: Path) -> set[str]:
    accessions = set()
    with open(str(INPUT_FASTA)) as in_handle:
        for seq_record in SeqIO.parse(in_handle, "fasta"):
            accessions.add(seq_record.id)

    return accessions


def get_GII_sequences(INPUT_METADATA: Path, accessions: set[str], genbank_mode = False) -> set[str]:
    gii_accessions: set[str] = set()

    if not genbank_mode:
        with open(INPUT_METADATA, "r", encoding="utf-8") as file:
            for line in file:
                meta_record: dict[str, Any] = json.loads(line)

                accession: str = meta_record["accession"]
                if accession in accessions:
                    if meta_record.get("virus", {}).get("taxId") == 122929 \
                        or meta_record.get("virus", {}).get("organismName") == "Norovirus GII":

                        gii_accessions.add(accession)

    elif genbank_mode:
        with INPUT_METADATA.open("r", encoding="utf-8") as in_handle:
            for seq_record in SeqIO.parse(in_handle, "genbank"):
                accession = seq_record.id

                if accession in accessions:
                    organism = seq_record.annotations.get("organism", "")
                    if organism == "Norovirus GII":
                        gii_accessions.add(accession)

    return gii_accessions


def get_sequence_completeness(INPUT_METADATA: Path) -> tuple[set[str], set[str]]:
    complete_genomes: set = set(); partial_genomes: set = set()
    
    with open(INPUT_METADATA, "r", encoding="utf-8") as file:
        for line in file:
            meta_record: dict[str, Any] = json.loads(line)

            accession: str = meta_record["accession"]
        
            match(meta_record["completeness"]):
                case("COMPLETE"):
                    complete_genomes.add(accession)
                case("PARTIAL"):
                    partial_genomes.add(accession)

    return (complete_genomes, partial_genomes)


def get_length_filtered_sequences(
    INPUT_FASTA: Path, 
    length_threshold: int
) -> set:
    cmd: list = ["seqkit", "seq", "-m", str(length_threshold), "-n", str(INPUT_FASTA)]
    
    length_filtered = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    sequences_above_threshold: set = {
        line.split()[0] 
        for line in length_filtered.stdout.splitlines() 
        if line.strip()
    }

    return sequences_above_threshold


def get_unique_sequences(INPUT_FASTA: Path) -> set[str]:
    cmd: list[str] = ["seqkit", "rmdup", "-i", "-n", str(INPUT_FASTA)]
    
    try:
        deduplicated_filtered = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as err:
        print("Command failed:")
        print("cmd", err.cmd)
        print("returncode", err.returncode)
        print("stdout", err.stdout)
        print("stderr", err.stderr)
        raise

    unique_sequences: set = {
        line[1:].split()[0] 
        for line in deduplicated_filtered.stdout.splitlines() 
        if line.startswith(">")
    }
 
    return unique_sequences


def get_annotated_sequences(INPUT_METADATA: Path) -> set:
    annotated_sequences: set = set()

    with open(INPUT_METADATA, "r", encoding="utf-8") as file:
        for line in file:
            meta_record = json.loads(line)

            accession: str = meta_record["accession"]

            if meta_record.get("isAnnotated") is True:
                annotated_sequences.add(accession)

    return annotated_sequences


def get_ambiguous_filtered_sequences(
    INPUT_FASTA: Path,
    ambiguous_threshold: int
) -> set[str]:
    passed_sequences: set[str] = set()
    allowed = set("ATCG")

    with open(str(INPUT_FASTA)) as in_handle:
        for seq_record in SeqIO.parse(in_handle, "fasta"):
            accession: str = str(seq_record.id)
            seq: str = str(seq_record.seq).upper()

            ambiguous_count: int = sum(1 for base in seq if base not in allowed)
            ambiguous_fraction: float = ambiguous_count / len(seq)

            if ambiguous_fraction < ambiguous_threshold:
                passed_sequences.add(accession)

    return passed_sequences


def extract_and_save_to_fasta(
    INPUT_FASTA: Path,
    accessions: set[str],
    OUTPUT_FASTA: Path,
) -> Path:
    with INPUT_FASTA.open("r", encoding="utf-8") as in_handle, \
         OUTPUT_FASTA.open("w", encoding="utf-8") as out_handle:

        for seq_record in SeqIO.parse(in_handle, "fasta"):
            accession = seq_record.id.split(":", 1)[0]

            if accession in accessions:
                SeqIO.write(seq_record, out_handle, "fasta")

    return OUTPUT_FASTA


def extract_and_save_to_genbank(
    INPUT_GENBANK: Path,
    accessions: set[str],
    OUTPUT_GENBANK: Path,
) -> Path:
    with INPUT_GENBANK.open("r", encoding="utf-8") as in_handle, \
         OUTPUT_GENBANK.open("w", encoding="utf-8") as out_handle:

        for seq_record in SeqIO.parse(in_handle, "genbank"):
            accession = seq_record.id

            if accession in accessions:
                SeqIO.write(seq_record, out_handle, "genbank")

    return OUTPUT_GENBANK


def get_features(feature_name: str) -> set[str]:
    match feature_name:
        case "vp1":
            return vp1_names
        case "rdrp":
            return rdrp_names | polyprotein_names
        case "junction":
            return vp1_names | rdrp_names | polyprotein_names
        case _:
            raise ValueError(f"Unknown feature_name: {feature_name}")


def get_feature_name(seq_record) -> str:
    feature = seq_record.description.replace(seq_record.id, "", 1).strip()
    feature = feature.split("[", 1)[0].strip()
    return feature


def extract_cds_feature_to_fasta(
    INPUT_CDS_FASTA: Path, 
    accessions: set[str],
    OUTPUT_FASTA: Path,
    feature_name: str,
) -> Path:
    # written: dict[str: str] = defaultdict()

    with INPUT_CDS_FASTA.open("r", encoding="utf-9") as in_handle, \
        OUTPUT_FASTA.open("w", encoding="utf-8") as out_handle:
        
        features = get_features(feature_name)

        for seq_record in SeqIO.parse(in_handle, "fasta"):
            accession = seq_record.id.split(":")[0]

            if accession not in accessions:
                continue

            feature = get_feature_name(seq_record)

            if feature not in features:
                continue

            SeqIO.write(seq_record, out_handle, "fasta")

    return OUTPUT_FASTA


def extract_junction_sequences_to_fasta(
    INPUT_FASTA: Path,
    REGION_CSV: Path,
    accessions: set[str],
    OUTPUT_FASTA: Path,
) -> Path:
    region_df = pd.read_csv(REGION_CSV).set_index("accession")

    with INPUT_FASTA.open("r", encoding="utf-8") as in_handle, \
        OUTPUT_FASTA.open("w", encoding="utf-8") as out_handle:

        for record in SeqIO.parse(in_handle, "fasta"):
            accession = record.id.split(":", 1)[0]

            if accession not in accessions:
                continue

            if accession not in region_df.index:
                continue

            row = region_df.loc[accession]

            rdrp_start = int(row["rdrp_start"])
            vp1_end = int(row["vp1_end"])

            start = rdrp_start - 1
            end = vp1_end

            junction_record = SeqRecord(
                record.seq[start:end],
                id=f"{accession}:{rdrp_start}-{vp1_end}",
                description="ORF1_ORF2_junction",
            )

            SeqIO.write(junction_record, out_handle, "fasta")
    
    return OUTPUT_FASTA


def typing_tool_intialise(INPUT_FASTA: Path, name: str, batch_size: int = 500):
    JOB_DIR: Path = INPUT_FASTA.parent / "temporary_web_crawler_data"
    JOB_DIR.mkdir(parents=True, exist_ok=True)

    JOB_STATE_TSV = JOB_DIR / f"{name}_job_ids.tsv"

    seq_count = sum(1 for _ in SeqIO.parse(str(INPUT_FASTA), "fasta"))
    n_batches: int = max(1, math.ceil(seq_count / batch_size))

    BATCHES_PATH: list[Path] = [
        JOB_DIR / f"{name}_batch_{i + 1}.fasta"
        for i in range(n_batches)
    ]

    with ExitStack() as stack:
        BATCH_HANDLES = [
            stack.enter_context(path.open("w", encoding="utf-8"))
            for path in BATCHES_PATH
        ]

        for i, seq_record in enumerate(SeqIO.parse(str(INPUT_FASTA), "fasta")):
            batch_index = i % n_batches
            SeqIO.write(seq_record, BATCH_HANDLES[batch_index], "fasta")

    URL = "https://mpf.rivm.nl/mpf/typingtool/norovirus/"
    job_ids = []

    print(BATCH_HANDLES)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        for FASTA in BATCHES_PATH:
            UPLOAD_FASTA = Path(FASTA).resolve()

            page.goto(URL, wait_until="domcontentloaded")

            file_input_selector = 'input[type="file"][name="data"]'
            page.set_input_files(file_input_selector, str(UPLOAD_FASTA)) 

            page.wait_for_function(
                """() => {
                    const el = document.querySelector('input[type="file"][name="data"]');
                    return el && el.files && el.files.length > 0;
                }""",
                timeout=10000,
            )

            page.wait_for_function(
                """() => {
                    const el = document.querySelector('span.error-text.error');
                    if (!el) return false;
                    const style = window.getComputedStyle(el);
                    return style.display === 'none';
                }""",
                timeout=30000,
            )
            
            # page.locator('button[id^="button_run_"]').click()
            page.get_by_role("button", name="Start!").click()

            page.wait_for_url("**/job/**", timeout=120000)

            job_url = page.url.rstrip("/")
            job_id = job_url.split("/")[-1]
            job_ids.append(job_id)

            print(f"Submitted {FASTA.name}: job {job_id}")

            time.sleep(2)

        browser.close()

    with open(JOB_STATE_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for job_id in job_ids:
            writer.writerow([job_id, "in_progress"])


def typing_tool_get_results(JOB_DIR: Path, name: str) -> Path | None:
    JOB_STATE_TSV = JOB_DIR / f"{name}_job_ids.tsv"

    job_ids_status = defaultdict()

    with open(JOB_STATE_TSV, newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            job_id, status = row[0], row[1]
            job_ids_status[job_id] = status

    URL = "https://mpf.rivm.nl/mpf/typingtool/norovirus/"

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        page.goto(URL, wait_until="domcontentloaded")

        for job_id, status in job_ids_status.items():
            if status == "completed":
                print(f"Job {job_id} has already been downloaded.")
                continue

            try:
                page.locator('input[type="text"][size="10"]').fill(str(job_id))
                page.get_by_role("button", name="Go!").click()

                csv_link = page.get_by_role("link", name="Table (CSV format)")
                csv_link.wait_for(timeout=30000)

                with page.expect_download(timeout=12000) as download_info:
                    csv_link.click()

                download = download_info.value
                output_path = JOB_DIR / f"{name}_job_{job_id}_table.csv"
                download.save_as(str(output_path))

                job_ids_status[job_id] = "completed"
                print(f"Job {job_id} sucessfully downloaded.")

                page.goto(URL, wait_until="domcontentloaded")

            except PlaywrightTimeoutError: 
                print(f"Job {job_id} not ready yet. Keeping as in_progress.")
                job_ids_status[job_id] = "in_progress"

                page.goto(URL, wait_until="domcontentloaded")
                continue

            time.sleep(2)

        browser.close()

    with open(JOB_STATE_TSV, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for job_id, status in job_ids_status.items():
            writer.writerow([job_id, status])

    download_complete = all(status == "completed" for status in job_ids_status.values())

    if download_complete == True:
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        COMBINED_CSV = JOB_DIR / f"{name}_typing_results_{timestamp}.csv"

        dfs = []

        for file in JOB_DIR.glob(f"*{name}_job_*_table.csv"):
            df = pd.read_csv(file)
            dfs.append(df)
    
        combinded = pd.concat(dfs, ignore_index=True)
        combinded.to_csv(COMBINED_CSV, index=False)

        return COMBINED_CSV
    
    return None


def clean_up_temporary_files(TARGET_DIR):
    shutil.rmtree(TARGET_DIR)


def get_genomic_info(
    OUTPUT_CSV: Path,
    typing_information: bool = True,
    TYPING_FILE: Path | None = None,
    region_information: bool = False,
    CDS_FASTA: Path | None = None,
    genbank_mode: bool = False,
    GENBANK_FILE: Path = None
) -> Path | None:
    if region_information and typing_information:
        region_dict = defaultdict(make_row_full_information)
    elif typing_information:
        region_dict = defaultdict(make_row_typing_information)
    elif region_information:
        region_dict = defaultdict(make_row_region_information)

    if region_information:
        if genbank_mode is False:
            if not CDS_FASTA:
                print(f"CDS_FASTA is required when region_information=True")

            for record in SeqIO.parse(CDS_FASTA, "fasta"):
                accession, coords = record.id.split(":")
                begin, end = coords.split("-")
                begin = int(begin)
                end = int(end)

                coding_region = record.description.replace(record.id, "", 1).strip()
                coding_region = coding_region.split("[", 1)[0].strip()
                
                if coding_region in rdrp_names:
                    region_dict[accession]["rdrp_start"] = begin
                    region_dict[accession]["rdrp_end"] = end
                    region_dict[accession]["rdrp_source"] = "annotated_rdrp"
                    continue

                if coding_region in polyprotein_names:
                    if region_dict[accession]["rdrp_start"] is None and end >= 5000:
                        region_dict[accession]["rdrp_start"] = max(begin, end - 2200)
                        region_dict[accession]["rdrp_end"] = end
                        region_dict[accession]["rdrp_source"] = "polyprotein_fallback"
                    continue

                if coding_region in vp1_names:
                    region_dict[accession]["vp1_start"] = begin
                    region_dict[accession]["vp1_end"] = end
                    continue

        elif genbank_mode is True:
            for seq_record in SeqIO.parse(GENBANK_FILE, "genbank"):
                accession = seq_record.id

                for feature in seq_record.features:
                    if feature.type not in ["CDS", "mat_peptide"]:
                        continue
                    
                    coding_region = feature.qualifiers.get("product", [""])[0]

                    if not coding_region:
                        continue

                    begin = int(feature.location.start) + 1
                    end = int(feature.location.end)

                    if coding_region in rdrp_names:
                        region_dict[accession]["rdrp_start"] = begin
                        region_dict[accession]["rdrp_end"] = end
                        region_dict[accession]["rdrp_source"] = "annotated_rdrp"
                        continue

                    if coding_region in polyprotein_names:
                        if region_dict[accession]["rdrp_start"] is None and end >= 5000:
                            region_dict[accession]["rdrp_start"] = begin
                            region_dict[accession]["rdrp_end"] = end
                            region_dict[accession]["rdrp_source"] = "polyprotein_fallback"
                        continue

                    if coding_region in vp1_names:
                        region_dict[accession]["vp1_start"] = begin
                        region_dict[accession]["vp1_end"] = end
                        continue

        for accession in region_dict.keys():
            rdrp_start = region_dict[accession]["rdrp_start"]
            rdrp_end = region_dict[accession]["rdrp_end"]
            vp1_start = region_dict[accession]["vp1_start"]
            vp1_end = region_dict[accession]["vp1_end"]

            if rdrp_end is not None and vp1_start is not None:
                overlap = rdrp_end - vp1_start + 1

                if overlap > 0:
                    region_dict[accession]["has_overlap"] = True
                    region_dict[accession]["overlap"] = overlap
                else:
                    region_dict[accession]["has_overlap"] = False
                    region_dict[accession]["overlap"] = 0
            else:            
                region_dict[accession]["has_overlap"] = False
                region_dict[accession]["overlap"] = 0

            
            if region_dict[accession]["rdrp_source"] == "polyprotein_fallback" \
                and (
                    rdrp_start is not None
                    and rdrp_end is not None
                    and rdrp_start < 3500
                    and rdrp_end > 5000
                ):
                region_dict[accession]["rdrp_start"] = 3500
    
    if typing_information:
        with open(TYPING_FILE, "r", encoding="utf-8", newline="") as csvfile:
            reader = csv.DictReader(csvfile)

            for row in reader:
                accession = row["name"]
                region_dict[accession]["BLAST_score"] = row["BLAST score"]
                region_dict[accession]["begin"] = row["begin"]
                region_dict[accession]["end"] = row["end"]
                region_dict[accession]["length"] = row["length"]

                p_type = row["polymerase type"].strip()
                if p_type == "":
                    p_type = None
                else:
                    p_type = p_type.split()[0]

                region_dict[accession]["p_type"] = p_type
                
                p_subtype = row["polymerase subtype"].strip()
                if p_subtype == "":
                    p_subtype = "None"
                region_dict[accession]["p_subtype"] = p_subtype

                genotype_subtype = row["capsid subtype"]
                if genotype_subtype == "":
                    genotype_subtype = "None"
                region_dict[accession]["genotype"] = row["capsid type"]
                region_dict[accession]["genotype_subtype"] = genotype_subtype

    if not region_dict:
        return None
    
    df = pd.DataFrame.from_dict(region_dict, orient="index")
    df.index.name = "accession"
    df = df.reset_index()

    df.to_csv(OUTPUT_CSV, index=False)

    return OUTPUT_CSV


# def export_excluded_sequences(OUTPUT_CSV: Path, **sets_dict: set[str]) -> None:
#     with open(OUTPUT_CSV, 'w', encoding="utf-8", newline='') as csvfile:
#         writer = csv.writer(csvfile)
#         writer.writerow([
#             "accession",
#             "not_gii",
#             "partial",
#             "under_length_threshold",
#             "duplicated",
#             "not_annotated",
#             "ambiguous"
#         ])  

#         for accession in sorted(sets_dict["excluded_accessions"]):
#             writer.writerow([
#                 accession,
#                 accession not in sets_dict["gii_sequences"],
#                 accession not in sets_dict["complete_sequences"],
#                 accession not in sets_dict["length_filtered_sequences"],
#                 accession not in sets_dict["unique_sequences"],
#                 accession not in sets_dict["annotated_sequences"],
#                 accession not in sets_dict["non_ambiguous_sequences"]
#             ])

def export_excluded_sequences(
    output_csv: Path,
    excluded_accessions: set[str],
    **criteria: set[str]
) -> None:
    with open(output_csv, "w", encoding="utf-8", newline="") as csvfile:
        writer = csv.writer(csvfile)

        writer.writerow(["accession", *criteria.keys()])

        for accession in sorted(excluded_accessions):
            writer.writerow([
                accession,
                *[accession not in passed_set for passed_set in criteria.values()]
            ])


def get_rdrp_sequences(REGION_CSV: Path, min_length: int = 100) -> set[str]:
    df = pd.read_csv(REGION_CSV)

    rdrp_length = df["rdrp_end"] - df["rdrp_start"] + 1

    rdrp_df = df[
        df["rdrp_start"].notna()
        & df["rdrp_end"].notna()
        & (rdrp_length >= min_length)
    ]

    return set(rdrp_df["accession"].astype(str))


def get_vp1_sequences(REGION_CSV: Path, min_length: int = 100) -> set[str]:
    df = pd.read_csv(REGION_CSV)

    vp1_length = df["vp1_end"] - df["vp1_start"] + 1

    vp1_df = df[
        df["vp1_start"].notna()
        & df["vp1_end"].notna()
        & (vp1_length >= min_length)
    ]

    return set(vp1_df["accession"].astype(str))


def get_junction_sequences(
    REGION_CSV: Path,
    min_rdrp_length: int = 100,
    min_vp1_length: int = 100
) -> set[str]:
    df = pd.read_csv(REGION_CSV)

    rdrp_length = df["rdrp_end"] - df["rdrp_start"] + 1
    vp1_length = df["vp1_end"] - df["vp1_start"] + 1

    junction_df = df[
        df["rdrp_start"].notna()
        & df["rdrp_end"].notna()
        & df["vp1_start"].notna()
        & df["vp1_end"].notna()
        & (rdrp_length >= min_rdrp_length)
        & (vp1_length >= min_vp1_length)
        & (df["has_overlap"] == True)
    ]

    return set(junction_df["accession"].astype(str))


def get_genbank_sequence_completeness(INPUT_GB: Path) -> tuple[set[str],set[str]]:
    partial_set = set()
    complete_set = set()

    for seq_record in SeqIO.parse(INPUT_GB, "genbank"):
        accession = seq_record.id

        has_partial_cds = False

        for feature in seq_record.features:
            if feature.type != "CDS":
                continue

            start = feature.location.start
            end = feature.location.end

            if isinstance(start, BeforePosition) or isinstance(end, AfterPosition):
                has_partial_cds = True
                break

        if has_partial_cds:
            partial_set.add(accession)
        else:
            complete_set.add(accession)

    return (complete_set, partial_set)


def merge_fasta(INPUT_FASTAS: Path, OUTPUT_FASTA: Path) -> Path:
    records = []

    for fasta in INPUT_FASTAS:
        records.extend(SeqIO.parse(fasta, "fasta"))

    SeqIO.write(records, OUTPUT_FASTA, "fasta")
    return OUTPUT_FASTA

def keep_sequences_with_rdrp_and_vp1(specification_csv: Path) -> set[str]:
    df = pd.read_csv(specification_csv, low_memory=False)

    required_cols = [
        "accession",
        "rdrp_start",
        "rdrp_end",
        "vp1_start",
        "vp1_end",
    ]

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"{specification_csv} is missing required columns: {missing_cols}"
        )

    filtered_df = df[
        df["accession"].notna()
        & df["rdrp_start"].notna()
        & df["rdrp_end"].notna()
        & df["vp1_start"].notna()
        & df["vp1_end"].notna()
    ].copy()

    return set(filtered_df["accession"].astype(str).str.strip())


def prepare_final_complete_specification(df: pd.DataFrame) -> pd.DataFrame:
    required_cols = [
        "accession",
        "rdrp_start",
        "rdrp_end",
        "vp1_start",
        "vp1_end",
    ]

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df = df.copy()
    df["accession"] = df["accession"].astype(str).str.strip()

    region_cols = [
        "rdrp_start",
        "rdrp_end",
        "vp1_start",
        "vp1_end",
    ]

    df["_region_score"] = df[region_cols].notna().sum(axis=1)

    df = (
        df
        .sort_values("_region_score")
        .drop_duplicates("accession", keep="last")
        .drop(columns="_region_score")
    )

    df = df[
        df["rdrp_start"].notna()
        & df["rdrp_end"].notna()
        & df["vp1_start"].notna()
        & df["vp1_end"].notna()
    ].copy()

    return df


def filter_csv_on_ptype_genotype(INPUT_CSV: Path, OUTPUT_CSV: Path) -> Path:
    df = pd.read_csv(INPUT_CSV)

    ptype_clean = df["p_type"].fillna("").astype(str).str.strip()
    genotype_clean = df["genotype"].fillna("").astype(str).str.strip()

    invalid_values = [""] #, "could", "could not assign"]

    has_ptype = ~ptype_clean.str.lower().isin(invalid_values)
    has_genotype = ~genotype_clean.str.lower().isin(invalid_values)

    filtered_df = df[has_ptype | has_genotype].copy()

    filtered_df.replace("Could", "Could not assign", inplace=True)
    
    filtered_df.to_csv(OUTPUT_CSV, index=False)

    return OUTPUT_CSV

def assert_accessions(INPUT_FASTA: Path, INPUT_CSV: Path) -> set[str]:
    fasta_accessions = set()

    with INPUT_FASTA.open("r", newline="") as in_handle:
        for seq_record in SeqIO.parse(in_handle, "fasta"):
            fasta_accessions.add(seq_record.id)
    
    csv_accessions = set()
    with INPUT_CSV.open("r", encoding="utf-8", newline="") as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            csv_accessions.add(row["accession"])
    
  
    return fasta_accessions ^ csv_accessions

def filter_csv_accessions(INPUT_CSV: Path, OUTPUT_CSV: Path, accessions: set[str]) -> Path:
    df = pd.read_csv(INPUT_CSV)

    normalized_accessions = {str(acc).strip() for acc in accessions}

    csv_accessions = df["accession"].fillna("").astype(str).str.strip()

    filtered_df = df[csv_accessions.isin(normalized_accessions)].copy()

    filtered_df.to_csv(OUTPUT_CSV, index=False)

    return OUTPUT_CSV
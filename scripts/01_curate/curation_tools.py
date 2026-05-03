import subprocess
import datetime
import pandas as pd
import shutil
import json
import math
import time
import csv

from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError
from typing import Any
from contextlib import ExitStack

from collections import defaultdict
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def get_GII_sequences(
    INPUT_FASTA: Path, 
    INPUT_METADATA: Path,
    OUTPUT_FASTA: Path,
) -> set[str]:
    gii_accessions: set[str] = set()

    with open(INPUT_METADATA, "r", encoding="utf-8") as file:
        for line in file:
            meta_record: dict[str, Any] = json.loads(line)

            accession: str = meta_record["accession"]
        
            if meta_record.get("virus", {}).get("taxId") == 122929 \
                or meta_record.get("virus", {}).get("organismName") == "Norovirus GII":

                gii_accessions.add(accession)

    with open(str(INPUT_FASTA)) as in_handle, \
         open(str(OUTPUT_FASTA), 'w', encoding="utf8") as out_handle:

        for seq_record in SeqIO.parse(in_handle, "fasta"):
            seq_record: SeqRecord

            accession: str = str(seq_record.id)

            if accession in gii_accessions:
                SeqIO.write(seq_record, out_handle, "fasta")

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
        # print(deduplicated_filtered.stdout)
        # print(deduplicated_filtered.stderr)
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
    output_name: str
) -> Path:
    OUTPUT_FASTA = Path(INPUT_FASTA.parent / f"{output_name}.fna")

    with open(INPUT_FASTA) as in_handle, \
        open(OUTPUT_FASTA, "w") as out_handle:

        for seq_record in SeqIO.parse(in_handle, "fasta"):
            accession = seq_record.id.split(":")[0]
            if accession in accessions:
                SeqIO.write(seq_record, out_handle, "fasta")

    return OUTPUT_FASTA


def typing_tool_intialise(INPUT_FASTA: Path):
    JOB_DIR: Path = INPUT_FASTA.parent / "temporary_web_crawler_data"
    JOB_DIR.mkdir(parents=True, exist_ok=True)

    JOB_STATE_TSV = JOB_DIR / "job_ids.tsv"

    seq_count = sum(1 for _ in SeqIO.parse(str(INPUT_FASTA), "fasta"))
    n_batches: int = max(1, math.ceil(seq_count / 500))

    BATCHES_PATH: list[Path] = [
        JOB_DIR / f"batch_{i + 1}.fasta"
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


def typing_tool_get_results(JOB_DIR) -> Path | None:
    JOB_STATE_TSV = JOB_DIR / "job_ids.tsv"

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
                csv_link.wait_for(timeout=3000)

                with page.expect_download(timeout=1200) as download_info:
                    csv_link.click()

                download = download_info.value
                output_path = JOB_DIR / f"job_{job_id}_table.csv"
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
        COMBINED_CSV = JOB_DIR.parent / f"typing_results_{timestamp}.csv"

        dfs = []

        for file in JOB_DIR.glob("*job_*_table.csv"):
            df = pd.read_csv(file)
            dfs.append(df)
    
        combinded = pd.concat(dfs, ignore_index=True)
        combinded.to_csv(COMBINED_CSV, index=False)

        return COMBINED_CSV
    
    return None


def clean_up_temporary_files(TARGET_DIR):
    shutil.rmtree(TARGET_DIR)


def get_genomic_region_info(
    TYPING_FILE: Path, 
    option: str = "typing", 
    CDS_FASTA: Path | None = None
) -> Path | None:
    def make_row_complete_region():
        return {
            "BLAST_score"     : None,
            "begin"           : None,
            "end"             : None,
            "length"          : None,
            "p_type"          : None,
            "p_subtype"       : None,
            "rdrp_start"      : None,
            "rdrp_end"        : None,
            "genotype"        : None,
            "genotype_subtype": None,
            "vp1_start"       : None,
            "vp1_end"         : None
        }
    
    def make_row_typing_region():
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

    if option == "typing":
        region_dict = defaultdict(make_row_typing_region)

    elif option == "complete":
        if not CDS_FASTA and not accessions:
            print(f"Provide supplementary files for 'complete' information: CDS_FASTA and accessions of intrest.")

        region_dict = defaultdict(make_row_complete_region)

        vp1_names = {
            "VP1",
            "capsid VP1",
            "viral protein 1",
            "capsid protein VP1",
            "major capsid protein",
            "major viral capsid protein",
            "major capsid protein VP1"
        }

        rdrp_names = {
            "RdRp",
            "RNA-dependent RNA polymerase"
        }

        polyprotein_names = {
            "polyprotein",
            "nonstructural polyprotein"
        }

        for record in SeqIO.parse(CDS_FASTA, "fasta"):
            accession, coords = record.id.split(":")

            begin, end = coords.split("-")

            coding_region = record.description.replace(record.id, "", 1).strip()
            coding_region = coding_region.split("[", 1)[0].strip()

            if coding_region in vp1_names:
                region_dict[accession]["vp1_start"] = int(begin)
                region_dict[accession]["vp1_end"] = int(end)
                continue
                
            if coding_region in rdrp_names:
                region_dict[accession]["rdrp_start"] = int(begin)
                region_dict[accession]["rdrp_end"] = int(end)
                continue

            if coding_region in polyprotein_names:
                if region_dict[accession]["rdrp_start"] is not None:
                    continue
                region_dict[accession]["rdrp_start"] = 3500
                region_dict[accession]["rdrp_end"] = int(end)
                continue
    else:
        print(f"Provide a valid option for information (default:'typing', 'complete').")
        return None
    
    with open(TYPING_FILE, "r", encoding="utf-8", newline="") as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            accession = row["name"]
            region_dict[accession]["BLAST_score"] = row["BLAST score"]
            region_dict[accession]["begin"] = row["begin"]
            region_dict[accession]["end"] = row["end"]
            region_dict[accession]["length"] = row["length"]

            p_subtype = row["polymerase subtype"]
            if p_subtype == "":
                p_subtype = "None"
            region_dict[accession]["p_type"] = row["polymerase type"].split()[0]
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

    OUTPUT_CSV = TYPING_FILE.parent / "specification.csv"
    df.to_csv(OUTPUT_CSV, index=False)

    return OUTPUT_CSV
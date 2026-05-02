import subprocess
import json
import yaml
import math

from playwright.sync_api import sync_playwright
from typing import Any

from collections import defaultdict
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def get_GII_sequences(INPUT_FASTA: Path, 
                      INPUT_METADATA: Path,
                     ) -> tuple[Path, set[str]]:
    
    gii_accessions: set[str] = set()

    with open(INPUT_METADATA, "r", encoding="utf-8") as file:
        for line in file:
            meta_record: dict[str, Any] = json.loads(line)

            accession: str = meta_record["accession"]
        
            if meta_record.get("virus", {}).get("taxId") != 122929 \
                or meta_record.get("virus", {}).get("organismName") == "Norovirus GII":

                gii_accessions.add(accession)

    OUTPUT_FASTA: Path = INPUT_FASTA.parent / "gii_sequences.fna"

    with open(str(INPUT_FASTA)) as in_handle, \
         open(str(OUTPUT_FASTA), 'w', encoding="utf8") as out_handle:
        

        for seq_record in SeqIO.parse(in_handle, "fasta"):
            seq_record: SeqRecord

            accession: str = str(seq_record.id)

            if accession in gii_accessions:
                SeqIO.write(seq_record, out_handle, "fasta")

    return (OUTPUT_FASTA, gii_accessions)


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


def get_length_filtered_sequences(INPUT_FASTA: Path, length_threshold: int) -> set:
    cmd: list = ["seqkit", "seq", "-m", length_threshold, "-n", str(INPUT_FASTA)]
    
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


def get_unique_sequences(INPUT_FASTA: Path) -> set:
    cmd: list = ["seqkit", "rmdup", "-ignore-case", "-n", str(INPUT_FASTA)]
    
    deduplicated_filtered = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    unique_sequences: set = {
        line.split()[0] 
        for line in deduplicated_filtered.stdout.splitlines() 
        if line.strip()
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


def get_ambigious_filtered_sequences(INPUT_FASTA: Path, ambigious_threshold: int) -> set:
    sequences_above_threshold: set = set()
    allowed = set("ATCG")


    with open(str(INPUT_FASTA)) as in_handle:
        for seq_record in SeqIO.parse(in_handle, "fasta"):
            seq_record: SeqRecord
            accession: str = str(seq_record.id)

            seq: str = str(seq_record.seq).upper()
            ambiguous_count: int = sum(1 for base in seq if base not in allowed)
            ambiguous_fraction: float = ambiguous_count / len(seq)

            if ambiguous_fraction > ambigious_threshold:
                sequences_above_threshold.add(accession)

    return sequences_above_threshold


def typing_tool_intialise(INPUT_FASTA):
    seq_count: int = 0
    accessions: set = set()

    with open(str(INPUT_FASTA)) as in_handle:
        for seq_record in SeqIO.parse(in_handle, "fasta"):
            accessions.add(seq_record.id) 
            seq_count += 1

    CURATE_DIR = Path(__file__).resolve().parent
    PROJECT_ROOT = (CURATE_DIR.parent).parent
    CONFIG_PATH = PROJECT_ROOT / "config" / "curation.yml"

    with open(CONFIG_PATH, "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)

    JOB_DIR: Path = PROJECT_ROOT / config["paths"]["temp_typing_dir"]
    JOB_DIR.mkdir(parents=True, exist_ok=True)

    n_batches: int = max(1, math.ceil(seq_count / 1000))

    BATCHES_PATH: list[Path] = [
        JOB_DIR / f"batch_{i + 1}.fasta" for i in range(n_batches)
    ]

    JOB_HANDLES = [
        open(batch_path, "w", encoding="utf-8") 
        for batch_path in BATCHES_PATH
    ]

    with open(INPUT_FASTA, "r", encoding="utf-8") as in_handle:
        for i, seq_record in enumerate(SeqIO.parse(in_handle, "fasta")):
            batch_index = i % n_batches
            SeqIO.write(seq_record, JOB_HANDLES[batch_index], "fasta")

    for handle in JOB_HANDLES:
        handle.close()


    URL = "https://mpf.rivm.nl/mpf/typingtool/norovirus/"

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        page.goto(URL, wait_until="domcontentloaded")

        file_input_selector = 'input[type="file"][name="data"]'
        page.set_input_files(file_input_selector, str(INPUT_FASTA)) 

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
        
        start_button_selector = 'button[id^="button_run_"]'
        page.locator(start_button_selector).click()

        page.wait_for_url("**/job/**", timeout=120000)

        job_url = page.url.rstrip("/")
        job_id = job_url.split("/")[-1]


def run_norovirus_typing_tool(INPUT_FASTA, OUTPUT_DIR):
    URL = "https://mpf.rivm.nl/mpf/typingtool/norovirus/"

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        page.goto(URL, wait_until="domcontentloaded")

        file_input_selector = 'input[type="file"][name="data"]'
        page.set_input_files(file_input_selector, str(INPUT_FASTA)) 

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
        
        start_button_selector = 'button[id^="button_run_"]'
        page.locator(start_button_selector).click()

        page.wait_for_url("**/job/**", timeout=120000)

        job_url = page.url.rstrip("/")
        job_id = job_url.split("/")[-1]

        csv_link_selector = 'a[id^="csv-table-download_"]'
        page.wait_for_selector(csv_link_selector, timeout=300000)

        with page.expect_download(timeout=120000) as download_info:
            page.locator(csv_link_selector).click()

        download = download_info.value
        output_path = OUTPUT_DIR / f"{INPUT_FASTA.stem}_job_{job_id}_table.csv"
        download.save_as(str(output_path))

        print(f"Downloaded CSV to: {output_path}")
        browser.close()
        
        return output_path


def get_typing_regions(INPUT_FASTA: Path) -> dict[str, Any] | None:
    
    def make_row():
        return {
            "BLAST_score": None,
            "length": None,
            "genotype": None,
            "p_type": None,
            "rdrp_start": None,
            "rdrp_end": None,
            "vp1_start": None,
            "vp1_end": None
        }

    region_dict = defaultdict(make_row)
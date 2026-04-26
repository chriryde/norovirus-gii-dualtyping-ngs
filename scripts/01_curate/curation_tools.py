import pandas as pd
import typing
import subprocess
import json
import yaml

from playwright.sync_api import sync_playwright

from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def get_GII_sequences(INPUT_FASTA: Path, 
                      INPUT_METADATA: Path,
                     ) -> tuple[Path, set[str]]:
    
    gii_accessions: set = set()

    with open(INPUT_METADATA, "r", encoding="utf-8") as file:
        for line in file:
            record: str = json.loads(line)

            accession: str = record["accession"]
        
            if record.get("virus", {}).get("taxId") != 122929 \
                or record.get("virus", {}).get("organismName") == "Norovirus GII":

                gii_accessions.add(accession)

    OUTPUT_FASTA: Path = INPUT_FASTA.parent() / "gii_sequences.fna"

    with open(str(INPUT_FASTA)) as in_handle, \
         open(str(OUTPUT_FASTA), 'w', encoding="utf8") as out_handle:
        
        for record in SeqIO.parse(in_handle, "fasta"):
            accession: str = record.id

            if accession in gii_accessions:
                SeqIO.write(record, out_handle, "fasta")

    return (OUTPUT_FASTA, gii_accessions)


def get_sequence_completeness(INPUT_METADATA: Path) -> tuple[set[str], set[str]]:
    complete_genomes: set = set(); partial_genomes: set = set()
    
    with open(INPUT_METADATA, "r", encoding="utf-8") as file:
        for line in file:
            record: str = json.loads(line)

            accession: str = record["accession"]
        
            match(record["completeness"]):
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
            record = json.loads(line)

            accession: str = record["accession"]

            if record.get("isAnnotated") is True:
                annotated_sequences.add(accession)

    return annotated_sequences


def get_ambigious_filtered_sequences(INPUT_FASTA: Path, ambigious_threshold: int) -> set:
    sequences_above_threshold: set = set()
    allowed = set("ATCG")

    with open(str(INPUT_FASTA)) as in_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            accession: str = record.id

            seq: str = (record.seq).upper()
            ambiguous_count: int = sum(1 for base in seq if base not in allowed)
            ambiguous_fraction: float = ambiguous_count / len(seq)

            if ambiguous_fraction > ambigious_threshold:
                sequences_above_threshold.add(accession)

    return sequences_above_threshold


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


def get_typing_regions(INPUT_FASTA: Path) -> dict[int, int, str, str, int, int, int ,int]:
    
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
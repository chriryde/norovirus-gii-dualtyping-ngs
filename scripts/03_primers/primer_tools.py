
import csv
import pandas as pd
import datetime
import math
import time
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from Bio import SeqIO


def extract_amplicons(
    LEFT_PRIMER_FILE: Path,
    left_primer: str,
    RIGHT_PRIMER_FILE: Path,
    right_primer: str,
    INPUT_FASTA: Path,
    OUTPUT_DIR: Path
) -> Path:
    
    left_primer = f"varVAMP_{left_primer}_LEFT"
    right_primer = f"varVAMP_{right_primer}_RIGHT"

    start = 0
    with open(LEFT_PRIMER_FILE, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        for row in tsv_reader:
            if left_primer == row[2]:
                start = int(row[5])
                break 
    
    stop = 0
    with open(RIGHT_PRIMER_FILE, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        for row in tsv_reader:
            if right_primer == row[2]:
                stop = int(row[6])
                break 

    assert(start != 0)
    assert(stop != 0)

    OUTPUT_FILE = OUTPUT_DIR / f"{left_primer}-{start}-{right_primer}-{stop}.fna"

    print("Extracting seq between:", start, stop)

    with open(OUTPUT_FILE, "w") as out:
        for record in SeqIO.parse(INPUT_FASTA, "fasta"):
            extracted_seq = record.seq[start:stop] ####

            assert(extracted_seq)
            out.write(f">{record.id}_region_{start}_{stop}\n")
            out.write(str(extracted_seq) + "\n")

    print(f"Amplicons saved in {OUTPUT_FILE}")

    return OUTPUT_FILE


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
                csv_link.wait_for(timeout=3000)

                with page.expect_download(timeout=1200) as download_info:
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


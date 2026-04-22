import pandas as pd
import subprocess
import json
import yaml

from collections import defaultdict
from playwright.sync_api import sync_playwright
from pathlib import Path
from Bio import SeqIO

def classify_sequences(INPUT_FASTA, OUTPUT_DIR):
    URL = "https://mpf.rivm.nl/mpf/typingtool/norovirus/"

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        page.goto(URL, wait_until="networkidle")

        page.set_input_files('input[type="file"][name="data"]', str(INPUT_FASTA))   
from Bio import SeqIO
import pandas as pd
from sys import argv
import correct_primer_positions as cpp
import yaml
from pathlib import Path
import csv

##change from original script; adding iupac lexion for degenerate bases:

#pick allowed number of mismathces here
allowed_n_mismathces = 2
max_l_dist = 2

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "typing.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

VP1_FASTA = PROJECT_ROOT/ Path(config['paths']['vp1_partial'])
OUTPUT_FASTA_VP1 = PROJECT_ROOT / Path(config["paths"]["output_fasta_vp1_amplicon_ny_kort"])

IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

#VP1_primer = 'CWGCWGTGAACGCRTTCCC'
VP1_primer = 'GCAAGCCCCTAATGGTGAGTTT'


def iupac_mismatches(primer, target):
    mismatches = 0

    for p, t in zip(primer.upper(), target.upper()):
        primer_allowed = IUPAC.get(p, {p})
        target_allowed = IUPAC.get(t, {t})

        if primer_allowed.isdisjoint(target_allowed):
            mismatches += 1

    return mismatches


def find_iupac_matches(primer, reference, max_mismatch=allowed_n_mismathces):
    matches = []
    primer_len = len(primer)
    reference = reference.upper()
    


    for i in range(0, len(reference) - primer_len + 1):
        window = reference[i:i + primer_len]
        mismatches = iupac_mismatches(primer, window)

        if mismatches <= max_mismatch:
            matches.append({
                "start": i,
                "end": i + primer_len,
                "mismatches": mismatches
                })

    return matches


def find_best_iupac_match(primer, reference):
    primer = primer.upper()
    reference = reference.upper()
    primer_len = len(primer)

    best = None

    for i in range(0, len(reference) - primer_len + 1):
        window = reference[i:i + primer_len]
        mismatches = iupac_mismatches(primer, window)

        if best is None or mismatches < best["mismatches"]:
            best = {
                "start": i,
                "end": i + primer_len,
                "mismatches": mismatches,
                "window": window
            }

    return best


def type_VP1(OUTPUT_FASTA_VP1, VP1_FASTA, primer):
    start = 0
    with open(OUTPUT_FASTA_VP1, "w") as out:
            for record in SeqIO.parse(VP1_FASTA, "fasta"):
                primer_coord = find_best_iupac_match(primer, record)['end']
                #if 330 < primer_coord < 350:
                extracted_seq = record.seq[start:primer_coord] ####
                    
                assert(extracted_seq)
                out.write(f">{record.id}_region_{start}_{primer_coord}\n")
                out.write(str(extracted_seq) + "\n")
            # else:
            #     print('sequence out of range')

type_VP1(OUTPUT_FASTA_VP1, VP1_FASTA, VP1_primer)
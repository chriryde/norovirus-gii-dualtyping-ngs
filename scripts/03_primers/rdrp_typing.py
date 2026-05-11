from Bio import SeqIO
import pandas as pd
from sys import argv
import os
import yaml
from pathlib import Path

PRIMER_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PRIMER_DIR.parent.parent

CONFIG_PATH = PROJECT_ROOT / "config" / "amplicon_typing.yml"

with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

INPUT_FASTA = PROJECT_ROOT/ Path(config['paths']['rdrp_partial'])
OUTPUT_FASTA = PROJECT_ROOT / Path(config["paths"]["output_fasta_rdrp"])


#pick allowed number of mismathces here
allowed_n_mismathces = 4


rdrp_primer = str("GGTGATGATGARATWGTSAGCAC")

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


def find_best_iupac_match(rdrp_primer, reference):
    rdrp_primer = rdrp_primer.upper()
    reference = reference.upper()
    rdrp_primer_len = len(rdrp_primer)

    best = None

    for i in range(0, len(reference) - rdrp_primer_len + 1):
        window = reference[i:i + rdrp_primer_len]
        mismatches = iupac_mismatches(rdrp_primer, window)

        if best is None or mismatches < best["mismatches"]:
            best = {
                "start": i,
                "end": i + rdrp_primer_len,
                "mismatches": mismatches,
                "window": window
            }

    return best

def extract_from_primer_to_end(input_fasta, rdrp_primer, output_fasta):

    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            genome_seq = str(record.seq).upper()

            matches = find_iupac_matches(
                rdrp_primer,
                genome_seq,
                max_mismatch=allowed_n_mismathces
            )

            if len(matches) == 1:
                start = matches[0]["start"]
                stop = len(genome_seq)
                extracted_seq = genome_seq[start:stop]

                out.write(f">{record.id}_region_{start}_to_{stop}\n")
                out.write(extracted_seq + "\n")

            elif len(matches) > 1:
                print(f"WARNING: {record.id}: primer found multiple times, skipping.")
                print(matches)

            else:
                print(f"WARNING: {record.id}: primer not found, skipping.")
                best = find_best_iupac_match(rdrp_primer, genome_seq)
                print(f"Best match: {best}")


if __name__ == "__main__":

    extract_from_primer_to_end(
        INPUT_FASTA,
        rdrp_primer,
        OUTPUT_FASTA
    )
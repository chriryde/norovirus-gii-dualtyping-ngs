import subprocess
from Bio import SeqIO
import pandas as pd


def cluster_sequences(input_file, output_file, cutoff = "0.98"):
    """
    Cluster sequences using cd-hit-est with a specified cutoff.
    input_file: Path to the input FASTA file containing sequences to be clustered.
    output_file: Path to the output FASTA file where clustered sequences will be saved.
    cutoff: Sequence identity threshold for clustering (default is 0.98).
    Returns the path to the output file containing clustered sequences.
    """
    try:
        print("Start running cd hit \n")
        with open(output_file, "w") as out_f:
            cd_hit_cmd = ["cd-hit-est", "-i", str(input_file), "-o", str(output_file), "-c", cutoff]
            cluster = subprocess.run(
                cd_hit_cmd,
                stderr=subprocess.PIPE,
                check=True
            )
    except subprocess.CalledProcessError as err:
        print("Command failed:")
        print("cmd", err.cmd)
        print("returncode", err.returncode)
        print("stdout", err.stdout)
        print("stderr", err.stderr)
        raise
    print(f"Process finished without errors\n")

    return (output_file)

#################################################


def get_accesions(input_file):
    """
    Function that extracts accession numbers from a FASTA file and return them as a set.
    input_file: Path to the input FASTA file containing sequences.
    Returns a set of accession numbers extracted from the FASTA file.
    """

    accession_set = set()
    with open(str(input_file)) as in_handle:

        for record in SeqIO.parse(in_handle, "fasta"):
            accession_set.add(record.id)

    return accession_set

#################################################


def get_sequences(input_file, output_file, accession_set):
    """Function that extracts sequences from a FASTA file based on a set of accession numbers and saves them to an output FASTA file.
    input_file: Path to the input FASTA file containing sequences.
    output_file: Path to the output FASTA file where extracted sequences will be saved.
    accession_set: A set of accession numbers to filter the sequences.
    Returns the path to the output file containing the extracted sequences."""
    
    with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            accession = record.id
            if accession in accession_set:
                SeqIO.write(record, out_handle, "fasta")

    return output_file

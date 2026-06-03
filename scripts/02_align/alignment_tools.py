import csv
from pathlib import Path
import yaml
import pandas as pd
import subprocess
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def alignment(input_file, output_file, log_file):
    """
    Function for running an MSA for a FASTA file using MAFFT.

    input_file: FASTA file containing sequences to be aligned
    output_file: aligned sequences in FASTA format
    log_file: file storing the log output from MAFFT
    """
    
    try:
        print("Start running mafft \n")
        with open(output_file, "w") as out_f, open(log_file, "w") as out_log:
            mafft_cmd = ["mafft", "--auto", str(input_file)]
            alignment = subprocess.run(
                mafft_cmd,
                stdout=out_f,
                # stderr=subprocess.PIPE,
                stderr=out_log,
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


    with open(log_file, "r") as in_handle:
        lines = in_handle.readlines()
        for index, line in enumerate(lines):
            if lines[index] == "Strategy:\n":
                print(lines[index])
                print(lines[index+1])
                print(lines[index+2])

    print(f'The following file has been aligned: {input_file.name}')
    print('--------\n')
    return (output_file, log_file)


########################################################

def get_region(input_metadata_csv, input_file, output_file, start, end, type):
    """"
    The function picks out a specific region from each sequence in a FASTA file, 
    using stored metadata. 

    input_metadata_csv: file containing metadata for each sequence
    input_file: FASTA file containing genomes
    output_file: new FASTA file containing only the desired region from each genome
    start: the starting position for the desired region
    end: the ending position for the desired region
    type: specifies if the regions is for the p_type or genotype
    """


    df = pd.read_csv(input_metadata_csv)


    seq_dict = {}

    for index, row in df.iterrows():

        if not math.isnan(row[start]) and not math.isnan(row[end]) and not row[type] == "Could" and not pd.isna(row[type]) and not row[type] == "Could not assign":
            seq_dict[row['accession']] = (int(row[start]), int(row[end]))



    with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            accession = record.id
            if accession in seq_dict.keys():
                seq = str(record.seq).upper()

                region = seq[seq_dict[accession][0]: seq_dict[accession][1]+1]
            
                new_record = SeqRecord(
                    Seq(region),
                    id = accession,
                    name = accession,
                    description = record.description
                )

                SeqIO.write(new_record, out_handle, "fasta")

    print(df[type].unique())
    return output_file

#########################################################


def get_genotype_or_ptype(input_file, input_metadata_csv, output_file, column, type):
    """
    Function for filtering sequences in a FASTA file based on their genotype or p-type, using stored metadata.

    input_file: FASTA file containing sequences to be filtered
    input_metadata_csv: file containing metadata for each sequence
    output_file: new FASTA file containing only the sequences of the desired genotype or p-type
    column: specifies if the filtering is for the p_type or genotype
    type: specifies the desired genotype or p-type
    """
    
    df = pd.read_csv(input_metadata_csv)

    seq_set = set()

    for index, row in df.iterrows():
       
        if row[column] == type:
                seq_set.add(row['accession'])
            
        

    with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:
            
        new_record = [record for record in SeqIO.parse(in_handle, "fasta") if record.id in seq_set]

        SeqIO.write(new_record, out_handle, "fasta")

    return output_file


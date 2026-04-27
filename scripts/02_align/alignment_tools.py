import csv
from pathlib import Path
import yaml
import pandas as pd
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


## funktion för alignment
def alignment(input_file, output_file, log_file):
    
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





## funktion för att plocka ut regioner
def get_region(input_metadata_csv, input_file, output_file, start, end):
    df = pd.read_csv(input_metadata_csv)


    seq_dict = {}

    # check that all names match the csv file
    for index, row in df.iterrows():
        ## update to the exact names in input_metadata_csv:
        ## update in pipe_align as well
        seq_dict[row['accession']] = (row[start], row[end])



    with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            accession = record.id
            seq = str(record.seq).upper()

            region = seq[seq_dict[accession][0]: seq_dict[accession][1]+1]
        
            new_record = SeqRecord(
                Seq(region),
                id = accession,
                name = accession,
                description = record.description
            )

            SeqIO.write(new_record, out_handle, "fasta")

    return output_file









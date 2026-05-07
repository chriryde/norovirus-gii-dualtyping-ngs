import subprocess
from Bio import SeqIO
import pandas as pd


def cluster_sequences(input_file, output_file, cutoff = "0.98"):
    try:
        print("Start running cd hit \n")
        with open(output_file, "w") as out_f:
            cd_hit_cmd = ["cd-hit-est", "-i", str(input_file), "-o", str(output_file), "-c", cutoff]
            cluster = subprocess.run(
                cd_hit_cmd,
                #stdout=out_f,
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


    


def get_accesions(input_file):
    accession_set = set()
    with open(str(input_file)) as in_handle:

        for record in SeqIO.parse(in_handle, "fasta"):
            accession_set.add(record.id)

    return accession_set

# if not math.isnan(row[start]) and not math.isnan(row[end]) and not row[type] == "Could" and not pd.isna(row[type]) and not row[type] == "Could not assign":

# function to get specific regions 
def get_sequences(input_file, output_file, accession_set):
    
    with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            accession = record.id
            if accession in accession_set:
                SeqIO.write(record, out_handle, "fasta")

    return output_file


# def remove_unsure_seq(input_file, output_file, metadata, accessions):
#     ## p_types kan vara "Could" eller nan
#     ## genotypes kan vara vara nan
#     df = pd.read_csv(metadata)
#     known_accessions = set()

#     # with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:

    
#     for index, row in df.iterrows():
#         if row["accession"] in accessions and not row["p_type"] == "Could" and not pd.isna(row["p_type"]) and not pd.isna(row["genotype"]):
#             known_accessions.add(row["accession"])

#     print(f'längd på known_accessions = {len(known_accessions)}')

#     return (get_sequences(input_file, output_file, known_accessions), known_accessions)


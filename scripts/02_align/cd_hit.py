import subprocess
from Bio import SeqIO

def cluster_sequences(input_file, output_file, cutoff = "0.98", wordsize = "10"):
    try:
        print("Start running cd hit \n")
        with open(output_file, "w") as out_f:
            cd_hit_cmd = ["cd-hit-est", "-i", str(input_file), "-o", str(output_file), "-c", cutoff, "-n", wordsize]
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



# function to get specific regions 
def get_sequences(input_file, output_file, accession_set):
    
    with open(str(input_file)) as in_handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            accession = record.id
            if accession in accession_set:
                SeqIO.write(record, out_handle, "fasta")

    return output_file

import pandas as pd
import cd_hit as cd
import alignment_tools as at
from Bio import SeqIO


def partial_accessions(complete_metadata_csv, partial_metadata_csv):
    df = pd.read_csv(complete_metadata_csv)


    genotype_dict = {}
    p_type_dict = {}
    could_not_assign = set()


    for index, row in df.iterrows():
       
        if row["p_type"] == "Could":
            could_not_assign.add(row['accession'])
        elif row["p_type"] in p_type_dict:
            p_type_dict[row["p_type"]] += 1
        else:
            p_type_dict[row["p_type"]] = 1
       


        if row["genotype"] == "Could not assign":
            could_not_assign.add(row['accession'])
        elif row["genotype"] in genotype_dict:
            genotype_dict[row["genotype"]] += 1
        else:
            genotype_dict[row["genotype"]] = 1


    df = pd.read_csv(partial_metadata_csv)

    rdrp_partial = set()
    vp1_partial = set()

   
    for index, row in df.iterrows():
        if row['p_type'] in p_type_dict.keys() and p_type_dict[row['p_type']] <= 15 and (row['rdrp_end']-row['rdrp_start']) > 1000:
                rdrp_partial.add(row['accession'])
       
        if row['genotype'] in p_type_dict.keys() and p_type_dict[row['genotype']] <= 15 and (row['vp1_end']-row['vp1_start']) > 1000:
                vp1_partial.add(row['accession'])

    return (p_type_dict, genotype_dict)


def m(partial_fasta, full_fasta, partial_metadata_csv, start, end):
    partial_output = at.get_region(partial_metadata_csv, partial_fasta, partial_output, start, end)


    with open("new_fasta", "w") as out_handle:
        records = list(SeqIO.parse(str(partial_output), "fasta")) + list(SeqIO.parse(str(full_fasta), "fasta"))
    
    SeqIO.write(records, out_handle, "fasta")

    return out_handle



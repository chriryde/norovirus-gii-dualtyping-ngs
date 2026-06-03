import pandas as pd
import cd_hit as cd
import alignment_tools as at
from Bio import SeqIO
import math


def partial_accessions(complete_metadata_csv, partial_metadata_csv):    
    """
    Takes in metadata_csv for complete and for partial. Then calculates how many of each genotype/p-type that occur in the complete metadata-set.
    Then iterates over the partial metadata. Returns set for rdrp and vp1 that contains accessions for the sequences that meet requirements.

    complete_metadata_csv: Path to the CSV file containing metadata for complete genomes.
    partial_metadata_csv: Path to the CSV file containing metadata for partial genomes.
    Returns a set of accessions for sequences that meet the specified requirements based on genotype/p-type"""
    

    df_full = pd.read_csv(complete_metadata_csv)

    genotype_dict = {}
    p_type_dict = {}
    could_not_assign = set()


    for index, row in df_full.iterrows():
       
        if row["p_type"] == "Could" or pd.isna(row["p_type"]):
            could_not_assign.add(row['accession'])
        elif row["p_type"] in p_type_dict:
            p_type_dict[row["p_type"]] += 1
        else:
            p_type_dict[row["p_type"]] = 1
       

        if row["genotype"] == "Could not assign" or pd.isna(row["genotype"]):
            could_not_assign.add(row['accession'])
        elif row["genotype"] in genotype_dict:
            genotype_dict[row["genotype"]] += 1
        else:
            genotype_dict[row["genotype"]] = 1


    
    df_partial = pd.read_csv(partial_metadata_csv)

    rdrp_partial = set()
    vp1_partial = set()
    common_dict = {}

   
    for index, row in df_partial.iterrows():
        if row['p_type'] in p_type_dict.keys() and p_type_dict[row['p_type']] <= 15 and (row['rdrp_end']-row['rdrp_start']) > 1000:
                
            rdrp_partial.add(row['accession'])

            if row["p_type"] in common_dict:
                    common_dict[row["p_type"]] += 1
            else:
                common_dict[row["p_type"]] = 1
                
        elif not row["p_type"] in p_type_dict.keys() and (row['rdrp_end']-row['rdrp_start']) > 1000:
            rdrp_partial.add(row['accession'])

            if row["p_type"] in common_dict:
                common_dict[row["p_type"]] += 1
            else:
                common_dict[row["p_type"]] = 1
       
        if row['genotype'] in genotype_dict.keys() and genotype_dict[row['genotype']] <= 15 and (row['vp1_end']-row['vp1_start']) > 1000:
                
            vp1_partial.add(row['accession'])


            if row["genotype"] in common_dict:
                common_dict[row["genotype"]] += 1
            else:
                common_dict[row["genotype"]] = 1
        elif not row["genotype"] in genotype_dict.keys() and (row['vp1_end']-row['vp1_start']) > 1000:
            vp1_partial.add(row['accession'])

            if row["genotype"] in common_dict:
                common_dict[row["genotype"]] += 1
            else:
                common_dict[row["genotype"]] = 1

            
    print('------------------------------------\n')
    return (rdrp_partial, vp1_partial)


#####################################

def merge_partial_regions(partial_fasta, full_fasta, partial_output):
    """Merges two FASTA files containing partial and full sequences, respectively, into a single FASTA file.
    partial_fasta: Path to the input FASTA file containing partial sequences.
    full_fasta: Path to the input FASTA file containing full sequences.
    partial_output: Path to the output FASTA file where merged sequences will be saved.
    Returns the path to the output file containing the merged sequences."""
    
    with open(partial_output, "w") as out_handle:
        records = list(SeqIO.parse(str(partial_fasta), "fasta")) + list(SeqIO.parse(str(full_fasta), "fasta"))
    
        SeqIO.write(records, out_handle, "fasta")

    return partial_output



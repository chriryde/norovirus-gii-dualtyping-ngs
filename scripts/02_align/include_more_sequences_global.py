import csv
from pathlib import Path
import yaml
import pandas as pd
import subprocess
import math


import cd_hit as cd
import alignment_tools as at

def global_accessions(complete_metadata_csv, partial_metadata_csv):

    """
    Takes in metadata_csv for complete and for partial. Then calculates how many of each genotype/p-type that occur in the complete metadata-set.
    Then iterates over the partial metadata. Returns set for rdrp and vp1 that contains accessions for the sequences that meet requirements.

    complete_metadata_csv: Path to the CSV file containing metadata for complete genomes.
    partial_metadata_csv: Path to the CSV file containing metadata for partial genomes.
    Returns a set of accessions for sequences that meet the specified requirements based on genotype/p-type
    """
    
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
    full_genomes_from_partial = set()
    common_dict = {}

   

    for index, row in df_partial.iterrows():
        if row['p_type'] in p_type_dict.keys() and p_type_dict[row['p_type']] <= 15 and row['length'] > 3500:
                
            full_genomes_from_partial.add(row['accession'])

            if row["p_type"] in common_dict:
                    common_dict[row["p_type"]] += 1
            else:
                common_dict[row["p_type"]] = 1
                
        elif not row["p_type"] in p_type_dict.keys() and row['length'] > 3500:
            full_genomes_from_partial.add(row['accession'])

            if row["p_type"] in common_dict:
                common_dict[row["p_type"]] += 1
            else:
                common_dict[row["p_type"]] = 1
                
        elif not row["p_type"] in p_type_dict.keys() and row['length'] > 3500:
            full_genomes_from_partial.add(row['accession'])

            if row["p_type"] in common_dict:
                common_dict[row["p_type"]] += 1
            else:
                common_dict[row["p_type"]] = 1

       
        if row['genotype'] in genotype_dict.keys() and genotype_dict[row['genotype']] <= 15 and row['length'] > 3500:              
            full_genomes_from_partial.add(row['accession'])


            if row["genotype"] in common_dict:
                common_dict[row["genotype"]] += 1
            else:
                common_dict[row["genotype"]] = 1
        elif not row["genotype"] in genotype_dict.keys() and row['length'] > 3500:
            full_genomes_from_partial.add(row['accession'])

            if row["genotype"] in common_dict:
                common_dict[row["genotype"]] += 1
            else:
                common_dict[row["genotype"]] = 1

            
    print(f'common_dict:\n{common_dict}\n---------------\n')
    print('------------------------------------\n')
    return (full_genomes_from_partial)

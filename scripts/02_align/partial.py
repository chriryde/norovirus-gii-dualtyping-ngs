import pandas as pd
import cd_hit as cd
import alignment_tools as at
from Bio import SeqIO
import math


def partial_accessions(complete_metadata_csv, partial_metadata_csv):
    """Tar in metadat_csv för complete och för partial. Beräknar sedan hur många av varje genotyp/p-typ som förekommer i commplete metadata-setet.
       Itererar sedan över partial metadata. Returnerar set för rdrp och vp1 som innheåller accessions för de skvenser som uppfyller krav."""
    
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



    print(f'alla förekommande genotyper i full: {genotype_dict}')
    print(f'alla förekommande p-typer i full: {p_type_dict}')
    
    df_partial = pd.read_csv(partial_metadata_csv)
    print(df_full["p_type"].unique())
    print(df_full["genotype"].unique())

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

            
    print(f'common_dict:\n{common_dict}\n---------------\n')

    print(f'antal accessions i båda seten = {len(rdrp_partial | vp1_partial)}')
    print(f'antal p-typer som borde läggas till = {len(rdrp_partial)}')
    print(f'antal genotyper som borde läggas till = {len(vp1_partial)}')
    print(f'antal sekvenser som borde tas med = {sum(common_dict.values())}')
    print('------------------------------------\n')
    return (rdrp_partial, vp1_partial)


def merge_partial_regions(partial_fasta, full_fasta, partial_output):

    with open(partial_output, "w") as out_handle:
        records = list(SeqIO.parse(str(partial_fasta), "fasta")) + list(SeqIO.parse(str(full_fasta), "fasta"))
    
        SeqIO.write(records, out_handle, "fasta")

    return partial_output



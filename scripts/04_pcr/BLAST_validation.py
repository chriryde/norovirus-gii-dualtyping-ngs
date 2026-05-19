import pandas as pd

def check_amplicons(input_csv):
    df = pd.read_csv(input_csv, sep="\t")
    
    genome_dict = {}
    primer_interaction_dict = {}

    for index, row in df.iterrows():
        
        if row['Distance'] == 'D' or row['Distance'] == 'D0':
            if row['Pairing and Strand direction ((+)/(-))'] not in primer_interaction_dict:
                primer_interaction_dict[row['Pairing and Strand direction ((+)/(-))']] = 1
            else:
                primer_interaction_dict[row['Pairing and Strand direction ((+)/(-))']] += 1
            


        
        elif row['Genome'] not in genome_dict:
            genome_dict[row['Genome']] = [1, 0]

            if row['Strand check'] == 'pass' and row['Amplification'] == 'exp' and row['Strand direction'] == '+-' and len(row['Distance']) > 1 and int(row['Distance'][1:]) > 600 and int(row['Distance'][1:]) < 900:
                genome_dict[row['Genome']][1] += 1

        else:
            genome_dict[row['Genome']][0] += 1

           
            if row['Strand check'] == 'pass' and row['Amplification'] == 'exp' and row['Strand direction'] == '+-' and len(row['Distance']) > 1 and int(row['Distance'][1:]) > 600 and int(row['Distance'][1:]) < 900: 
                genome_dict[row['Genome']][1] += 1

    return (genome_dict, primer_interaction_dict)


def passed_amplicons(data, output_file):
    rows = []

    for sequence, values in data.items():
        percent = (values[1]/values[0])*100

        rows.append([sequence, values[0], values[1], percent])

    df = pd.DataFrame(rows, columns=["Accession", "Total occurance", "Amplicons", "Percent of passed amplicons"])

    df.to_csv(output_file, index = False)

    return df


def passed_types(amplicon_csv, specification_csv, output_csv):

    df_amplicon = pd.read_csv(amplicon_csv)
    passed_accessions = []
    for index, row in df_amplicon.iterrows():
        if float(row['Percent of passed amplicons']) > 0:
            passed_accessions.append(row['Accession'])

    
    df_spec = pd.read_csv(specification_csv)
    passed_rows = []

    for index, row in df_spec.iterrows():
        if row['accession'] in passed_accessions:
            passed_rows.append(row)
    
    df_new = pd.DataFrame(passed_rows, columns=df_spec.columns)
    df_new.to_csv(output_csv, index = False)

    return df_new







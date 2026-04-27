
from Bio import AlignIO
from Bio.SeqUtils import seq1
import math

# Reads a MSA in FASTA format. 
alignment = AlignIO.read("alignment.fasta", "fasta")


def shannon_entropy_column(column):
    """
    Calculate the Shannon Entropy for a single column in the alignment.
    Lower entropy indicates higher conservation.
    """
    bases = {}
    total_count = len(column)

    # Count the occurrence of each nucleotide (or gap)
    for nt in column:
        if nt in bases:
            bases[nt] += 1
        else:
            bases[nt] = 1

    # Calculate entropy
    entropy = 0.0
    for count in bases.values():
        frequency = count / total_count
        entropy -= frequency * math.log(frequency, 2) # Log base 2

    return entropy

# Calculate entropy for every column in the alignment
def shannon_entropy_aligment(alignment_complete):
    all_columns = alignment_complete.get_alignment_length()
    return [shannon_entropy_column(alignment_complete[:, i]) for i in range(all_columns)]
    

def conserved_positions(alignment, threshold):
    entropy_scores = shannon_entropy_aligment(alignment)
    # Find the indices of conserved columns
    return [i for i, entropy in enumerate(entropy_scores) if entropy <= threshold]



def find_conserved_regions(alignment, threshold):
    region_min = 16
    conserved_col = conserved_positions(alignment, threshold)
    
    regions = []
    if not conserved_col:
        return regions
    start = conserved_col[0]

    for i in range(1, len(conserved_col)):
        # om inte sammanhängande
        if conserved_col[i] != conserved_col[i-1] + 1:
            end = conserved_col[i-1]
            
            if end - start + 1 >= region_min:
                regions.append((start, end))
            
            start = conserved_col[i]

    # sista regionen
    end = conserved_col[-1]
    if end - start + 1 >= region_min:
        regions.append((start, end))

    return regions

print(find_conserved_regions(alignment, 0.2))




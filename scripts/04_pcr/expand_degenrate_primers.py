from itertools import product
from Bio import SeqIO
from Bio.Data import IUPACData
import argparse

def expand_degenerate_sequence(seq: str):
    seq = seq.upper().replace("U", "T")

    iupac = IUPACData.ambiguous_dna_values

    bases_per_position = []

    for base in seq:
        if base not in iupac:
            raise ValueError(f"Invalid base '{base}' in sequence: {seq}")
        bases_per_position.append(iupac[base])

    for combo in product(*bases_per_position):
        yield "".join(combo)

def main():
    parser = argparse.ArgumentParser(
        description="Expand degenerate primers in FASTA-format"
    )

    parser.add_argument("-i", "--input", required=True, help="Input FASTA with degenrate primers")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA with expanded primers")
    parser.add_argument(
        "--max-variants",
        type=int,
        default=10000,
        help="Stop if a primer create more variants than 10000"
    )

    args = parser.parse_args()

    with open(args.output, "w") as out:
        for record in SeqIO.parse(args.input, "fasta"):
            seq = str(record.seq).upper().replace("U", "T")

            variant_count = 1
            for base in seq:
                variant_count *= len(IUPACData.ambiguous_dna_values[base])

            if variant_count > args.max_variants:
                raise ValueError(
                    f"{record.id} ger {variant_count} varianter, vilket är över gränsen "
                    f"{args.max_variants}. Höj --max-variants om du verkligen vill expandera."
                )

            for i, expanded_seq in enumerate(expand_degenerate_sequence(seq), start=1):
                out.write(f">{record.id}_var{i}\n{expanded_seq}\n")

            print(f"{record.id}: {variant_count} varianter")

if __name__ == "__main__":
    main()
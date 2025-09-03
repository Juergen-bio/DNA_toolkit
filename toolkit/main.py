
import argparse
from .dna_toolkit import *
from .structures import *
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description='Analyze DNA sequences for ORFs.')
    parser.add_argument('-f', '--file', help='Path to a FASTA file containing DNA sequences.')
    parser.add_argument('-s', '--sequence', help='Raw DNA sequence string to analyze.')

    args = parser.parse_args()

    dna_sequences = []

    if args.file:
        try:
            with open(args.file, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    dna_sequences.append(str(record.seq).upper())
        except FileNotFoundError:
            print(f"Error: The file '{args.file}' was not found.")
            return
    elif args.sequence:
        dna_sequences.append(args.sequence)
    else:
        print("Error: No file or sequence provided. Use -f or -s.")
        return

    if not dna_sequences:
        print("No DNA sequences found to process.")
        return

    concatenated_dna = ''.join(dna_sequences)

    # Get all six reading frames (3 forward, 3 reverse)
    all_frames = get_all_reading_frames(concatenated_dna, codon_table)

    # Find ORFs in each reading frame
    all_proteins = []
    for frame in all_frames:
        all_proteins.extend(reading_frames(frame))

    # Print the list of found ORFs
    for protein in all_proteins:
        print(protein)


if __name__ == "__main__":
    main()
'''

python main.py --file my_data.fasta
python main.py --sequence ATGCATGCGCGA
'''
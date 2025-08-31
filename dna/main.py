
from DNAToolkit import *
from structures import *
from Bio import SeqIO

dna_file = r'path'
dna_sequences = []

with open(dna_file, 'r') as handle:
  for record in SeqIO.parse(handle, 'fasta'):
    dna_sequences.append(str(record.seq).upper())

# because the function is expecting a single str and not multiple
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
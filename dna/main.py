
from DNAToolkit import *
from structures import *
from Bio import SeqIO

folder_path = r'path'
dna = []

with open(folder_path, 'r') as handle:
  for record in SeqIO.parse(handle, 'fasta'):
    dna.append(str(record.seq))

# because the function is expecting a single str and not multiple
concatenated_dna = ''.join(dna)

aa_seq = (translation(concatenated_dna))
print(reading_frames(aa_seq))


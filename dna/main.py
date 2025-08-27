from DNAToolkit import *
from structure import *
from Bio import SeqIO

folder_path = r'C:\Users\hp\PycharmProjects\pythonProject\project\maize_files\MON810.fasta'
dna = []

with open(folder_path, 'r') as handle:
  for record in SeqIO.parse(handle, 'fasta'):
    dna.append(str(record.seq))

# because the function is expecting a single str and not multiple
concatenated_dna = ''.join(dna)

aa_seq = (translation(concatenated_dna))
print(reading_frames(aa_seq))


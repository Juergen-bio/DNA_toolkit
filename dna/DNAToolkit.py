from structures import *

#check to validate it is a dna string regardless of the source.
def validate(dna_seq):
    temp = dna_seq.upper()
    for nuc in temp:
        if nuc not in nucleotide:
            return False
    return temp

def count(dna_seq):
    number = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for nuc in dna_seq:
        if nuc in number:
            number[nuc] += 1
    return number

def transcription(seq):
    return seq.replace('T', 'U')

def reversecompliment(seq):
    return ''.join([dna_reversecompliment[nuc] for nuc in seq])[::-1]

def gc_content(seq):
    return (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100


def translation(seq):
    aa_seq = ''
    for i in range (0, len(seq) - len(seq)% 3, 3):
        codon = seq[i:i + 3]
        aminoacid = codon_table[codon]
        if codon not in codon_table:
            aminoacid = 'X'
        else:
            aa_seq += aminoacid
    return aa_seq

def reading_frames(aa_seq):
    current_protein = []
    protein = []

    for aa in aa_seq:
        if aa == 'M': # starts a string of orf proteins
            if current_protein:
                protein.append(current_protein)
            current_protein = [aa]
        elif aa == '_':
            if current_protein:
                protein.append(current_protein)
                current_protein = []
        else:
            if current_protein:
                current_protein.append(aa)

    if current_protein: # adds any remaining proteins
        protein.append(current_protein)
    return protein

def orf_proteins(seq):
    reading = []
    reading.append(translation(seq, 0))
    reading.append(translation(seq, 1))
    reading.append(translation(seq, 2))
    reading.append(translation(reversecompliment(seq, 0)))
    reading.append(translation(reversecompliment(seq, 1)))
    reading.append(translation(reversecompliment(seq, 2)))


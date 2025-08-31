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
        if nuc.upper() in number:
            number[nuc.upper()] += 1
    return number

def transcription(seq):
    return seq.upper().replace('T', 'U')

def reverse_complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_seq = seq[::-1].upper()
    complement_seq = ''

    for nuc in reverse_seq:
        if nuc in complement_dict:
            complement_seq += complement_dict[nuc]
        else:
            complement_seq += 'N' #'N' as a placeholder for non-nucleotides

    return complement_seq

def gc_content(seq):
    return (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100


def translation(seq, codon_table):
    aa_seq = ''
    for i in range(0, len(seq) - len(seq)% 3, 3):
        codon = seq[i:i + 3].upper()

        if codon in codon_table:
            aa_seq += codon_table[codon]

        else:
            aa_seq += 'X'
    return aa_seq

def reading_frames(aa_seq):
    current_protein = []
    proteins = []

    for aa in aa_seq:
        if aa == 'M':  # Starts a new potential ORF
            if current_protein and current_protein[-1] != '_':
                proteins.append(''.join(current_protein))
            current_protein = [aa]
        elif aa == '_':  # Terminates a protein
            if current_protein:
                current_protein.append(aa)
                proteins.append(''.join(current_protein))
                current_protein = []
        else:
            if current_protein:
                current_protein.append(aa)

    if current_protein:
        proteins.append(''.join(current_protein))
    return proteins

def orf_proteins(seq, codon_table):
    reading = []
    reading.append(translation(seq, codon_table))
    reading.append(translation(seq, codon_table))
    reading.append(translation(seq, codon_table))
    reading.append(translation(reverse_complement(seq, codon_table)))
    reading.append(translation(reverse_complement(seq, codon_table)))
    reading.append(translation(reverse_complement(seq, codon_table)))
    return frames

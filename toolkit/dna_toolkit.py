from .structures import *
import re

# check to validate it is a dna string regardless of the source.
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
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].upper()

        if codon in codon_table:
            aa_seq += codon_table[codon]

        else:
            aa_seq += 'X'
    return aa_seq

def reading_frames(aa_seq):
    orf_list = []
    # Regular expression to find sequences starting with 'M' and ending with a stop codon '_'
    # A positive lookahead ensures the stop codon is part of the match, but not consumed, allowing
    # for chained ORFs in the same frame.
    pattern = re.compile(r'M(?!_)[^_]*_')

    for match in re.finditer(pattern, aa_seq):
        orf = match.group()
        # Ensure ORF is not empty and is not just a stop codon
        if len(orf) > 1:
            orf_list.append(orf)
    return orf_list

def get_all_reading_frames(seq, codon_table):
    frames = []
    # Forward frames
    for i in range(3):
        frames.append(translation(seq[i:], codon_table))
    # Reverse frames
    rev_comp = reverse_complement(seq)
    for i in range(3):
        frames.append(translation(rev_comp[i:], codon_table))
    return frames
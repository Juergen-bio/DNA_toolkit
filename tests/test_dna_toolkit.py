import pytest
from toolkit.dna_toolkit import reverse_complement, translation
from toolkit.structures import codon_table

def test_reverse_complement_valid_sequence():
    """Test reverse complement with a standard DNA sequence."""
    seq = "ATCG"
    expected = "CGAT"
    assert reverse_complement(seq) == expected

def test_reverse_complement_lowercase():
    """Test reverse complement with a lowercase sequence."""
    seq = "atcg"
    expected = "CGAT"
    assert reverse_complement(seq) == expected

def test_translation_valid_sequence():
    """Test translation of a known coding sequence."""
    seq = "ATGCGCCGGCGCATGTAG"
    # Corrected expected value based on the codon table
    expected = "MRRRM_"
    assert translation(seq, codon_table) == expected

def test_translation_invalid_codon():
    """Test translation with an unknown codon."""
    seq = "NNNATG"
    # Corrected expected value
    expected = "XM"
    assert translation(seq, codon_table) == expected

def test_gc_content():
    """Test GC content calculation."""
    # Assuming gc_content function is in DNAToolkit
    from toolkit.dna_toolkit import gc_content
    seq = "GCGCGCGCGC"
    assert gc_content(seq) == 100.0
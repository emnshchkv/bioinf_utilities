from typing import Union


def is_dna(sequence: str) -> bool:
    """
    Determines if sequence is DNA.

    Arguments:
    sequence: str

    Returns:
    True: if sequence is DNA
    False: if sequence is NOT DNA
    """
    dna_alphabet = set("ATGCatgc")
    unique_nucleotides = set(sequence)
    if unique_nucleotides <= dna_alphabet:
        return True
    return False

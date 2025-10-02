from typing import Union


def is_rna(sequence: str) -> bool:
    """
    Determines if sequence is RNA.

    Arguments:
    sequence: str

    Returns:
    True: if sequence is RNA
    False: if sequence is NOT RNA
    """
    rna_alphabet = set("AUGCaugc")
    unique_nucleotides = set(sequence)
    if unique_nucleotides <= rna_alphabet:
        return True
    return False

from . import is_rna, is_dna


def is_nucleic_acid(sequence: str) -> bool:
    """
    Verifies if the sequence is RNA or DNA.

    Arguments:
    sequence: str

    Returns:
    True: if sequence is a nucleic acid
    False: if sequence is NOT a nucleic acid
    """
    if is_rna.is_rna(sequence):
        return True
    elif is_dna.is_dna(sequence):
        return True
    return False

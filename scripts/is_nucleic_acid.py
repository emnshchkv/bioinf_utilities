import which_nucleic_acid


def is_nucleic_acid(sequence: str) -> bool:
    """
    Verifies if the sequence is RNA or DNA.

    Arguments:
    sequence: str

    Returns:
    True: if sequence is a nucleic acid
    False: if sequence is NOT a nucleic acid
    """
    if which_nucleic_acid(sequence) is None:
        return False
    else:
        return True

from . import is_nucleic_acid
from . import complement
from typing import Union


def reverse_complement(sequence: str) -> Union[str, None]:
    """
    Converts a sequence into its reverse-complement counterpart.

    Arguments:
    sequence: str

    Returns:
    str: sequence reverse-complemet counterpart
    None: if sequence is NOT a nucleic acid
    """
    if is_nucleic_acid.is_nucleic_acid(sequence):
        comp_seq = complement.complement(sequence)
        return comp_seq[::-1]
    else:
        return None

from . import is_nucleic_acid, complement, reverse
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
        return complement.complement(reverse.reverse(sequence))

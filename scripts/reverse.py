from . import is_nucleic_acid
from typing import Union


def reverse(sequence: str) -> Union[str, None]:
    """
    Converts a sequence into its reverse counterpart.

    Arguments:
    sequence: str

    Returns:
    str: reversed sequence
    None: if sequence is NOT a nucleic acid
    """
    if is_nucleic_acid.is_nucleic_acid(sequence):
        return sequence[::-1]

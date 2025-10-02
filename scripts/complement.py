from . import is_rna, is_dna
from typing import Union


def complement(sequence: str) -> Union[str, None]:
    """
    Converts a sequence into its complement counterpart.

    Arguments:
    sequence: str

    Returns:
    str: sequence complemet counterpart
    None: if sequence is NOT a nucleic acid

    """
    if is_dna.is_dna(sequence):
        comp_dict = {
            "A": "T",
            "a": "t",
            "T": "A",
            "t": "a",
            "G": "C",
            "g": "c",
            "C": "G",
            "c": "g",
        }
    elif is_rna.is_rna(sequence):
        comp_dict = {
            "A": "U",
            "a": "u",
            "U": "A",
            "u": "a",
            "G": "C",
            "g": "c",
            "C": "G",
            "c": "g",
        }
    else:
        return None
    countepart = ""
    for nucleotide in sequence:
        countepart += comp_dict[nucleotide]
    return countepart

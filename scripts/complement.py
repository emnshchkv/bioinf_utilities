import which_nucleic_acid
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
    if which_nucleic_acid.which_nucleic_acid(sequence) == "DNA":
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
        countepart = ""
        for nucleotide in sequence:
            countepart += comp_dict[nucleotide]
        return countepart
    elif which_nucleic_acid.which_nucleic_acid(sequence) == "RNA":
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
        countepart = ""
        for nucleotide in sequence:
            countepart += comp_dict[nucleotide]
        return countepart
    else:
        return None

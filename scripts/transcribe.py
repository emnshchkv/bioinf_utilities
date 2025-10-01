from . import is_nucleic_acid
from typing import Union


def transcribe(
    sequence: str,
    transcription_dict: dict = {
        "A": "U",
        "a": "u",
        "T": "A",
        "t": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
    },
) -> Union[str, None]:
    """
    Generates a trancsript for a DNA sequence.

    Arguments:
    sequence: str
    transcription_dict: dict

    Returns:
    str: for a DNA sequence
    None: for a RNA sequence
    """
    if is_nucleic_acid.is_nucleic_acid(sequence):
        mrna = ""
        for nucleotide in sequence:
            mrna += transcription_dict[nucleotide]
        return mrna
    else:
        return None

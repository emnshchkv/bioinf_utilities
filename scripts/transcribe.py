from . import is_dna
from typing import Union


def transcribe(sequence: str) -> Union[str, None]:
    """
    Generates a mRNA-trancsript for a DNA sequence.

    Arguments:
    sequence: str

    Returns:
    str: for a DNA sequence
    None: for a RNA sequence
    """
    transcription_dict: dict = {
        "A": "U",
        "a": "u",
        "T": "A",
        "t": "a",
        "G": "C",
        "g": "c",
        "C": "G",
        "c": "g",
    }
    if is_dna.is_dna(sequence):
        return "".join([transcription_dict[n] for n in sequence])

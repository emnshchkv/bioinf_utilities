from typing import Union


def is_dna(sequence: str) -> bool:
    """
    Determines if sequence is DNA.

    Arguments:
    sequence: str

    Returns:
    True: if sequence is DNA
    False: if sequence is NOT DNA
    """
    dna_alphabet = set("ATGCatgc")
    unique_nucleotides = set(sequence)
    return unique_nucleotides <= dna_alphabet


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
    return unique_nucleotides <= rna_alphabet


def is_nucleic_acid(sequence: str) -> bool:
    """
    Verifies if the sequence is RNA or DNA.

    Arguments:
    sequence: str

    Returns:
    True: if sequence is a nucleic acid
    False: if sequence is NOT a nucleic acid
    """
    return is_rna(sequence) or is_dna(sequence)


def reverse(sequence: str) -> Union[str, None]:
    """
    Converts a sequence into its reverse counterpart.

    Arguments:
    sequence: str

    Returns:
    str: reversed sequence
    None: if sequence is NOT a nucleic acid
    """
    if is_nucleic_acid(sequence):
        return sequence[::-1]


def complement(sequence: str) -> Union[str, None]:
    """
    Converts a sequence into its complement counterpart.

    Arguments:
    sequence: str

    Returns:
    str: sequence complemet counterpart
    None: if sequence is NOT a nucleic acid

    """
    if is_dna(sequence):
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
    elif is_rna(sequence):
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


def reverse_complement(sequence: str) -> Union[str, None]:
    """
    Converts a sequence into its reverse-complement counterpart.

    Arguments:
    sequence: str

    Returns:
    str: sequence reverse-complemet counterpart
    None: if sequence is NOT a nucleic acid
    """
    if is_nucleic_acid(sequence):
        return complement(reverse(sequence))


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
    if is_dna(sequence):
        return "".join([transcription_dict[n] for n in sequence])

def which_nucleic_acid(
    sequence: str, dna_alphabet: set =set("ATGCatgc"), rna_alphabet: set=set("AUGCaugc")
) -> str:
    """
    Determines the type of nucleic acid for a sequence.
    
    Arguments:
    sequence: str
    dna_alphabet: set
    rna_alphabet: set

    Returns:
    str: if sequence is a nucleic acid
    None: if sequence is NOT a nucleic acid
    """
    unique_nucleotides = set(sequence)
    if unique_nucleotides <= dna_alphabet:
        return "DNA"
    elif unique_nucleotides <= rna_alphabet:
        return "RNA"
    else:
        return None
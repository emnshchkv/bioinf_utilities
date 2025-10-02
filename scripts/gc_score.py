def gc_score(sequence: str) -> float:
    """
    Calculates GC-content score of a sequence.

    Arguments:
    sequence: str

    Returns:
    float: GC-content
    """
    sequence = sequence.lower()
    return (sequence.count("g") + sequence.count("c")) * 100 / len(sequence)

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


def quality_score(quality_line: str) -> float:
    """
    Calculates quatily score for a pair (sequence, quality).

    Arguments:
    quality_line: str

    Returns:
    float: average sequence quality
    """
    encoding = dict()
    for x in range(0, 41):
        encoding.update({chr(x + 33): int(x)})
    total_score = 0
    for sym in quality_line:
        total_score += encoding[sym]
    return total_score / len(quality_line)


def is_suitable_fastq(
    sequence: str,
    quality_line: str,
    gc_lower: int | float,
    gc_upper: int | float,
    len_lower: int | float,
    len_upper: int | float,
    quality_threshold: int | float,
) -> bool:
    """
    Determines if sequence should be kept based on filtering criteria.

    Arguments:
    sequence: str
    quality: str
    gc_lower: int | float
    gc_upper: int | float
    len_lower: int
    len_upper: int
    quality_thresh: int | float

    Returns:
    True: if sequence passes all filters
    False: if sequence DOES DOT pass all filters
    """
    gc_content = gc_score(sequence)
    sequence_quality = quality_score(quality_line)
    sequence_length = len(sequence)

    return (
        gc_lower <= gc_content <= gc_upper
        and len_lower <= sequence_length <= len_upper
        and sequence_quality >= quality_threshold
    )

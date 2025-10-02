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

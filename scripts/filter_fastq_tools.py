import os


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


def read_fastq(*input_path: str) -> dict:
    """
    Reads *.fastq and converts into 'identifier: (sequence, quality)' dictionary.

    Arguments:
    *input_path: absolute path

    Returns:
    dict[str, tuple[str, str]]
    """
    fastq_sequences = dict()
    with open(os.path.join(*input_path), "r") as fastq:
        while True:
            identifier = fastq.readline().strip()
            if not identifier:
                break
            sequence = fastq.readline().strip()
            _ = fastq.readline()
            quality = fastq.readline().strip()
            fastq_sequences[identifier] = (sequence, quality)
    return fastq_sequences


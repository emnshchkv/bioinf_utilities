from scripts import dna_rna_tools, filter_fastq_tools, assistant_ulitities
from typing import Union


def run_dna_rna_tools(*funnel: str) -> Union[str, list[str | None], bool]:
    """
    Processing DNA or RNA sequnce.

    Arguments:
    funnel : tuple[str]

    Returns:
    str: if single sequence was given
    list: if several sequences were given

    Raises error:
    if tool is unknown
    if no agruments are given
    """
    if len(funnel) < 2:
        raise ValueError(
            "At least one sequence and the name of the operation must be given."
        )
    *sequences, tool_name = funnel
    tools = {
        "is_rna": dna_rna_tools.is_rna,
        "is_dna": dna_rna_tools.is_dna,
        "is_nucleic_acid": dna_rna_tools.is_nucleic_acid,
        "transcribe": dna_rna_tools.transcribe,
        "reverse": dna_rna_tools.reverse,
        "complement": dna_rna_tools.complement,
        "reverse_complement": dna_rna_tools.reverse_complement,
    }
    if tool_name not in tools.keys():
        raise KeyError(
            f"Unknown tool: {tool_name}. Available tools: {list(tools.keys())}"
        )
    current_tool = tools[tool_name]
    result = [current_tool(seq) for seq in sequences]
    return result[0] if len(result) == 1 else result


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: int | tuple[int | float, int | float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    """
    Filtering fastq with manual settings.

    Agruments:
    seqs: dict
    gc_bounds: int | tuple
    length_bounds: int | tuple
    quality_threshold: int

    Returns:
    dict: filtered fastq
    """
    filtered = dict()
    (gc_lower_bound, gc_upper_bound), (len_lower_bound, len_upper_bound) = (
        assistant_ulitities.make_interval(gc_bounds),
        assistant_ulitities.make_interval(length_bounds),
    )
    for name, data in seqs.items():
        if (
            (gc_lower_bound <= filter_fastq_tools.gc_score(data[0]) <= gc_upper_bound)
            and (quality_threshold <= filter_fastq_tools.quality_score(data[1]))
            and (len_lower_bound <= len(data[0]) <= len_upper_bound)
        ):
            filtered[name] = data
    return filtered

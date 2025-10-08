import os
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
    fastq_path: str,
    gc_bounds: int | float | tuple[int | float, int | float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int | float = 0,
) -> None:
    """
    Filters *.fastq with manual settings and save it to 'output.fastq'.

    Agruments:
    fastq_path: absolute path to *.fastq
    gc_bounds: int | tuple
    length_bounds: int | tuple
    quality_threshold: int

    Returns:
    None
    """
    (gc_lower_bound, gc_upper_bound), (len_lower_bound, len_upper_bound) = (
        assistant_ulitities.make_interval(gc_bounds),
        assistant_ulitities.make_interval(length_bounds),
    )
    input_path = os.path.join(fastq_path)
    output_path = assistant_ulitities.make_output_path(fastq_path, "output.fastq")

    with open(input_path, "r") as input_fastq, open(output_path, "a") as output_fastq:
        while True:
            identifier = input_fastq.readline().strip()
            if not identifier:
                break
            sequence = input_fastq.readline().strip()
            _ = input_fastq.readline()
            quality = input_fastq.readline().strip()
            if filter_fastq_tools.is_suitable_fastq(
                sequence,
                quality,
                gc_lower_bound,
                gc_upper_bound,
                len_lower_bound,
                len_upper_bound,
                quality_threshold,
            ):
                output_fastq.write(
                    identifier + "\n" + sequence + "\n" + "+" + "\n" + quality
                )

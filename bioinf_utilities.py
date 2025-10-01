from scripts import (
    which_nucleic_acid,
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement,
)


def run_dna_rna_tools(*funnel: str):
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
    sequences = funnel[:-1]
    tool_name = funnel[-1]
    tools = {
        "which_nucleic_acid": which_nucleic_acid.which_nucleic_acid,
        "is_nucleic_acid": is_nucleic_acid.is_nucleic_acid,
        "transcribe": transcribe.transcribe,
        "reverse": reverse.reverse,
        "complement": complement.complement,
        "reverse_complement": reverse_complement.reverse_complement,
    }
    if tool_name not in tools.keys():
        raise KeyError(
            f"Unknown tool: {tool_name}. Available tools: {list(tools.keys())}"
        )
    current_tool = tools[tool_name]
    result = [current_tool(seq) for seq in sequences]
    return result[0] if len(result) == 1 else result

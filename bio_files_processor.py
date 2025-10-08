import os
from scripts import assistant_ulitities


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str | None = None
) -> None:
    """
    Converts *.fasta with multi-line sequences into *.fasta with one-line sequences.

    Arguments:
    input_fasta: absolute path to input FASTA
    output_fasta: absolute path to output FASTA

    Returns:
    None
    """
    if output_fasta is None:
        output_path = assistant_ulitities.make_output_path(input_fasta, "output.fasta")
    input_path = os.path.join(input_fasta)
    with open(input_path, "r") as input_fasta, open(output_path, "a") as output_fasta:
        current_sequence = []
        name = None
        for line in input_fasta:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    output_fasta.write(name + "\n" + "".join(current_sequence) + "\n")
                name = line
                current_sequence = []
            else:
                current_sequence.append(line)
        if name is not None:
            output_fasta.write(name + "\n" + "".join(current_sequence) + "\n")

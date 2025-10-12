import os, re
from scripts import assistant_ulitities


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str | None = None
) -> None:
    """
    Converts *.fasta with multi-line sequences into *.fasta with one-line sequences.

    Arguments:
    input_fasta: absolute path to input FASTA-file
    output_fasta: absolute path to output FASTA-file

    Returns:
    None
    """
    if output_fasta is None:
        output_path = assistant_ulitities.make_output_path(input_fasta, "output.fasta")
    else:
        output_path = os.path.join(output_fasta)
    input_path = os.path.join(input_fasta)
    with open(input_path, "r") as input_fasta_file, open(
        output_path, "a"
    ) as output_fasta_file:
        current_sequence = []
        name = None
        for line in input_fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    output_fasta_file.write(
                        name + "\n" + "".join(current_sequence) + "\n"
                    )
                name = line
                current_sequence = []
            else:
                current_sequence.append(line)
        if name is not None:
            output_fasta_file.write(name + "\n" + "".join(current_sequence) + "\n")


def parse_blast_output(input_file: str, output_file: str | None = None) -> None:
    """
    Parses the BLAST_output_file.txt and extracts the first protein from the Description column for each QUERY.

    Arguments:
    input_file: absolute path to BLAST_output_file.txt
    output_file: absolute path to the extracted proteins

    Returns:
    None
    """
    if output_file is None:
        output_path = assistant_ulitities.make_output_path(
            input_file, "output_parse_blast.txt"
        )
    else:
        output_path = os.path.join(output_file)
    input_path = os.path.join(input_file)
    with open(input_path, "r") as blast_output, open(
        output_path, "a"
    ) as parsed_blast_output:
        query = ""
        placed_in_section, table_start = False, False
        for line in blast_output:
            line = line.strip()
            if line.startswith("Query #"):
                query = line
                continue
            if "Sequences producing significant alignments:" in line:
                placed_in_section = True
                continue
            if "Description" in line:
                table_start = True
                continue
            if query and placed_in_section and table_start:
                columns = re.split(r"\s{2,}|\.{3,}", line)
                parsed_blast_output.write(query + "\n" + columns[0].strip() + "\n")
                query = ""
                placed_in_section, table_start = False, False
            continue

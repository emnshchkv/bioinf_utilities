from typing import Union
import os, re


def make_interval(
    obj: int | float | tuple[int | float, int | float],
) -> Union[tuple[int, int], tuple[int, float], tuple[float, float]]:
    """
    Defines interval boundaries.

    Arguments:
    obj: int - upper boundary, tuple[int, int] - lower and upper boundaries

    Returns:
    range: numbers within interval
    """
    if isinstance(obj, int) or isinstance(obj, float):
        return (0, obj)
    else:
        return (obj[0], obj[1])


def make_output_path(input_path: str, file_name: str) -> str:
    """
    Creates a directory for output file based on the path of the input file.

    Arguments:
    input_path: str
    file_name: str

    Returns:
    str: full path to the output file with the specified name
    """
    output_dir = os.path.dirname(input_path)
    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, file_name)


def extract_gene_and_protein(input_file: str) -> list[tuple[str]]:
    """
    Extracts gene names and their corresponding protein sequences from GBK-file.

    Arguments:
    inout_file: absolute path to the GBK-file

    Returns:
    A list of tuples where each tuple contains two strings:
        - gene_name (str): name of the gene extracted from '/gene=' qualifier
        - protein_sequence (str): complete amino acid sequence of the protein encoded
    """
    gene_protein = []
    sequence = ""
    gene_name = ""
    reading_translation = False
    translation_lines = []
    for line in input_file:
        line = line.strip()
        if "/gene=" in line:
            matches = re.findall(r'gene="([^"]+)"', line)
            if matches:
                gene_name = matches[0]
        elif "/translation=" in line:
            reading_translation = True
            translation_lines = []
            match = re.search(r'translation="([^"]*)', line)
            if match:
                translation_part = match.group(1)
                translation_lines.append(translation_part)
                if line.strip().endswith('"'):
                    reading_translation = False
                    sequence = "".join(translation_lines)
                    if gene_name and sequence:
                        gene_protein.append((gene_name, sequence))
                        sequence = ""
                        gene_name = ""
        elif reading_translation:
            if '"' in line:
                end_part = line.split('"')[0]
                translation_lines.append(end_part)
                reading_translation = False
                sequence = "".join(translation_lines)
                if gene_name and sequence:
                    gene_protein.append((gene_name, sequence))
                    sequence = ""
                    gene_name = ""
            else:
                translation_lines.append(line)
        elif line and r"\s{2,}(\w+)" in line:
            if gene_name and sequence and not reading_translation:
                gene_protein.append((gene_name, sequence))
                sequence = ""
                gene_name = ""
    return gene_protein

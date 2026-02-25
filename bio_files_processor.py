import re
from Bio import SeqIO


def convert_multiline_fasta_to_oneline(input_path: str, output_path: str) -> None:
    """
    Convertы a FASTA file with multiline sequences into a FASTA file
    where each sequence is written on a single line.

    Args:
        input_path (str): Path to the input FASTA file (multiline format).
        output_path (str): Path to the output FASTA file (one-line format).

    Returns:
        None: The converted FASTA file is written to output_path.
    """
    with open(input_path, "r") as multiline, open(output_path, "w") as oneline:
        records = SeqIO.parse(multiline, "fasta")
        SeqIO.write(records, oneline, "fasta")


def parse_blast_output(input_path: str, output_path: str) -> None:
    """
    Parses the BLAST_output_file.txt and extracts the first protein from the Description column for each QUERY.

    Args:
        input_file: absolute path to BLAST_output_file.txt
        output_path: absolute path to the extracted proteins file

    Returns:
        None
    """

    with open(input_path, "r") as blast_output, open(
        output_path, "w"
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


# AI wrote the following function
def select_genes_from_gbk_to_fasta(
    input_path: str,
    genes: tuple[str],
    output_path: str,
    n_before: int = 1,
    n_after: int = 1,
) -> None:
    """
    Extract amino acid sequences of CDS features neighboring specified genes
    from a GenBank file and write them to a FASTA file.

    The genes of interest themselves are excluded from the output.
    Neighboring genes are determined based on their order among CDS features
    within each record of the GenBank file.

    Args:
        input_path (str): Path to the input GenBank (.gbk/.gbff) file.
        genes (tuple[str]): Gene names (as specified in the /gene qualifier) whose neighboring CDS features should be extracted.
        output_path (str): Path to the output FASTA file.
        n_before (int) default=1: Number of CDS features to include upstream (before) each gene of interest.
        n_after (int) default=1: Number of CDS features to include downstream (after) each gene of interest.

    Returns:
        None: Writes the selected amino acid sequences to the specified FASTA file.
    """
    records = list(SeqIO.parse(input_path, "genbank"))
    cds_features = []
    for record in records:
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", [None])[0]
                translation = feature.qualifiers.get("translation", [None])[0]
                if gene_name and translation:
                    cds_features.append((gene_name, translation))

    interest_indexes = []

    for i, (gene_name, _) in enumerate(cds_features):
        if gene_name in genes:
            interest_indexes.append(i)

    neighbours_indexes = set()

    for index in interest_indexes:
        start = max(0, index - n_before)
        stop = min(len(cds_features), index + n_after + 1)

        for i in range(start, stop):
            if i != index:
                neighbours_indexes.add(i)

    with open(output_path, "w") as output_fasta_file:
        for i in sorted(neighbours_indexes):
            gene_name, translation = cds_features[i]
            output_fasta_file.write(f">{gene_name}\n{translation}\n")

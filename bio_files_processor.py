import os, re
from scripts import assistant_ulitities
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


def parse_blast_output(input_file: str, output_file: str | None = None) -> None:
    """
    Parses the BLAST_output_file.txt and extracts the first protein from the Description column for each QUERY.

    Arguments:
    input_file: absolute path to BLAST_output_file.txt
    output_file: absolute path to the extracted proteins file

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


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    n_before: int = 1,
    n_after: int = 1,
    *genes: str | tuple[str],
    output_fasta: str = "output.fasta"
):
    """
    Extracts amino acid sequences adjacent to the genes of interest from the GBK-file. The genes of interest are not displayed.

    Arguments:
    input_gbk: absolute path to GBK-file
    genes: genes of interest which neighbors are being searched
    n_before, n_after: the number of genes before and after the gene of interest
    output_fasta: name of the output FASTA-file

    Returns:
    None
    """
    output_path = assistant_ulitities.make_output_path(input_gbk, output_fasta)
    input_path = os.path.join(input_gbk)
    with open(input_path, "r") as input_gbk_file:
        gene_protein = assistant_ulitities.extract_gene_and_protein(input_gbk_file)
    interest_indexes = []
    for gene in genes:
        for pair_gene_protein in gene_protein:
            if gene == pair_gene_protein[0]:
                interest_indexes.append(gene_protein.index(pair_gene_protein))
    neighbours_indexes = []
    for index in interest_indexes:
        start = index - n_before
        stop = index + n_after + 1
        if index - n_before < 0:
            start = 0
        elif index + n_after + 1 > len(gene_protein):
            stop = len(gene_protein)
        for i in range(start, stop):
            neighbours_indexes.append(i)
        neighbours_indexes.remove(index)
    with open(output_path, "a") as output_fasta_file:
        for index in neighbours_indexes:
            output_fasta_file.write(
                ">" + gene_protein[index][0] + "\n" + gene_protein[index][1] + "\n"
            )

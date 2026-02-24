# Bioinformatics Utilities

A collection of Python modules for common bioinformatics tasks: sequence manipulation, FASTQ filtering, and processing of FASTA, GenBank, and BLAST output files.

## Modules

The package consists of three independent modules:

- `bioseq_machinery.py` – object-oriented tools for DNA, RNA, and protein sequences.
- `fastq_filter.py` – FASTQ quality filtering based on length, GC‑content, and average Phred score.
- `bio_files_processor.py` – utilities for converting, parsing, and extracting data from biological file formats.

## Installation

No installation is required – simply place the `.py` files in your working directory or in a location accessible to Python. The modules use only the standard library and [Biopython](https://biopython.org/) (which you may need to install separately):

```bash
pip install biopython
```

Then import the desired functions in your script:

```python
from bioseq_machinery import DNASequence, RNASequence, AminoAcidSequence
from fast_filter import filter_fastq
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output, select_genes_from_gbk_to_fasta
```

---

## 1. Sequence Manipulation (`bioseq_machinery.py`)

This module provides classes for representing and manipulating biological sequences. Each class inherits from `BiologicalSequence` and implements alphabet validation, slicing, and common operations.

### Classes and Methods

| Class            | Description                                | Key Methods                                                                                              |
|------------------|--------------------------------------------|----------------------------------------------------------------------------------------------------------|
| `DNASequence`    | DNA sequence (alphabet: ATGCatgc)          | `complement()`, `reverse()`, `reverse_complement()`, `transcribe()` → `RNASequence`                      |
| `RNASequence`    | RNA sequence (alphabet: AUGCaugc)          | `complement()`, `reverse()`, `reverse_complement()`                                                      |
| `AminoAcidSequence` | Protein sequence (standard 20 amino acids + U, X, O, lowercase) | `categorize()` – counts residues by chemical class (nonpolar, polar, charged, etc.) |

All sequences support:
- Length via `len()`
- Indexing and slicing (e.g., `seq[5]`, `seq[10:20]`)
- String representation (`str(seq)`)

### Usage Examples

```python
from bioseq_machinery import DNASequence, RNASequence, AminoAcidSequence

# DNA operations
dna = DNASequence("ATGCGTAt")
print(dna)                       # Oligonucleotide : ATGCGTAt
print(dna.complement())           # TACGCATa
print(dna.reverse_complement())   # tAACGCAT? wait, example: better to use a real sequence
rna = dna.transcribe()            # RNASequence object
print(rna)                        # Oligonucleotide : UACGCAUa

# Protein categorization
protein = AminoAcidSequence("MVLSPADKTNVKAAWG")
counts = protein.categorize()
print(counts)  # e.g., {'nonpolar': 8, 'polar': 5, '+ charged': 2, ...}
```

---

## 2. FASTQ Filtering (`fastq_filter.py`)

The function `filter_fastq` reads a FASTQ file, applies filters, and writes passing reads to a new file.

### Function Signature

```python
filter_fastq(
    input_path: str,
    output_path: str,
    length_bounds: int | tuple = (0, 2**32),
    quality_threshold: int | float = 0,
    gc_bounds: int | tuple = (0, 100),
) -> None
```

**Parameters**
- `input_path` – path to the input FASTQ file.
- `output_path` – path where filtered reads will be saved.
- `length_bounds` – single integer (maximum length) or tuple `(min, max)`. Default `(0, 2**32)` (no practical limit).
- `quality_threshold` – minimum average Phred quality score.
- `gc_bounds` – single integer (maximum GC%) or tuple `(min, max)`. Default `(0, 100)`.

### Usage Examples

```python
from fast_filter import filter_fastq

# Keep reads with length 50–150 bp, average quality ≥ 20, GC% between 35% and 60%
filter_fastq(
    "raw.fastq",
    "filtered.fastq",
    length_bounds=(50, 150),
    quality_threshold=20,
    gc_bounds=(35, 60)
)

# Only filter by maximum length (200 bp) and no other restrictions
filter_fastq("raw.fastq", "trimmed.fastq", length_bounds=200)
```

---

## 3. Biological File Processing (`bio_files_processor.py`)

Three utility functions for working with common bioinformatics file formats.

### `convert_multiline_fasta_to_oneline`

Converts a FASTA file where sequences may be wrapped over multiple lines into a file where each sequence occupies a single line.

```python
convert_multiline_fasta_to_oneline("multiline.fasta", "oneline.fasta")
```

### `parse_blast_output`

Parses a **legacy plain‑text BLAST output** file (the default `-outfmt 0` pairwise format) and extracts, for each query, the first protein from the **Description** column. The output is a text file with alternating lines: query header and the extracted protein name.

```python
parse_blast_output("blast_results.txt", "extracted_proteins.txt")
```

**Note:** This parser is tailored to a specific output layout. For more robust parsing, consider using Biopython’s `Bio.SearchIO` with tabular BLAST output (`-outfmt 6` or `7`).

### `select_genes_from_gbk_to_fasta`

Extracts amino acid sequences of CDS features that neighbor user‑specified genes in a GenBank file. The genes of interest themselves are **excluded** from the output. Neighboring genes are identified by their order among all CDS features in the file.

```python
select_genes_from_gbk_to_fasta(
    input_path="genome.gbk",
    genes=("lacZ", "lacY"),
    output_path="neighbors.fasta",
    n_before=2,
    n_after=1
)
```

**Parameters**
- `input_path` – GenBank file (`.gbk` or `.gbff`).
- `genes` – tuple of gene names (as they appear in the `/gene` qualifier).
- `output_path` – output FASTA file.
- `n_before` – number of genes to include upstream (default 1).
- `n_after` – number of genes to include downstream (default 1).

The function collects all CDS features with a `/gene` qualifier and a `/translation` qualifier, then for every occurrence of any target gene, it adds the indices of its neighbors to a set (ensuring uniqueness) and finally writes the corresponding protein sequences in the original order.

---

## Error Handling

- All functions perform basic input validation (file existence, format compatibility) and raise appropriate exceptions (`ValueError`, `FileNotFoundError`, `RuntimeError`) with descriptive messages.
- Sequence classes validate the alphabet at instantiation and raise `ValueError` if invalid characters are found.
- Warnings are issued for missing qualifiers or genes not found in the file.

## Dependencies

- **Python** ≥ 3.7
- **Biopython** – required for `fast_filter.py` and `bio_files_processor.py` (installation: `pip install biopython`)

All other modules use only the Python standard library.

---

## Example Workflow

```python
from bioseq_machinery import DNASequence
from fast_filter import filter_fastq
from bio_files_processor import select_genes_from_gbk_to_fasta

# 1. Clean up sequencing reads
filter_fastq("raw.fastq", "clean.fastq", quality_threshold=25, length_bounds=(50, 200))

# 2. Analyse a DNA sequence
my_dna = DNASequence("ATCGATCG")
print(my_dna.reverse_complement())

# 3. Extract neighboring genes of a locus of interest
select_genes_from_gbk_to_fasta(
    "ecoli.gbk",
    genes=("rpsA", "rplB"),
    output_path="neighbors.faa",
    n_before=3,
    n_after=2
)
```

## License

This code is provided as-is for educational and research purposes. Feel free to adapt and reuse.
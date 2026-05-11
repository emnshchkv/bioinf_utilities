# Bioinformatics Utilities

A collection of Python modules for common bioinformatics tasks: sequence manipulation, FASTQ filtering, and processing of FASTA, GenBank, and BLAST output files.

## Modules

The package consists of 3 independent modules + test suite:

- `bioseq_machinery.py` – object-oriented tools for DNA, RNA, and protein sequences.
- `fastq_filter.py` – FASTQ quality filtering based on length, GC‑content, and average Phred-score.
- `bio_files_processor.py` – utilities for converting, parsing, and extracting data from biological file formats.
- `fastq_filter_tests.py` – pytest test suite for `fastq_filter.py`.

## Installation

No installation is required – simply place the `.py` files in your working directory or in a location accessible to Python. The modules depend on [biopython](https://biopython.org/) and [loguru](https://github.com/Delgan/loguru), which you can install with:

```bash
pip install biopython loguru
```

For running the tests, also install pytest:

```bash
pip install pytest
```

Then import the desired functions in your script:

```python
from bioseq_machinery import DNASequence, RNASequence, AminoAcidSequence
from fastq_filter import filter_fastq
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
print(dna)                        # Oligonucleotide : ATGCGTAt
print(dna.complement())           # TACGCATa
print(dna.reverse_complement())   # aTACGCAT
rna = dna.transcribe()            # RNASequence object
print(rna)                        # Oligonucleotide : UACGCAUa

# Protein categorization
protein = AminoAcidSequence("MVLSPADKTNVKAAWG")
counts = protein.categorize()
print(counts)  # e.g., {'nonpolar': 8, 'polar': 5, '+ charged': 2, ...}
```

---

## 2. FASTQ Filtering (`fastq_filter.py`)

The function `filter_fastq` reads a FASTQ file, applies filters, and writes passing reads to a new file. Progress and results are logged via `loguru` to both the console and a log file.

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
- `quality_threshold` – minimum average Phred quality score. Default `0`.
- `gc_bounds` – single integer (maximum GC%) or tuple `(min, max)`. Default `(0, 100)`.

### Usage Examples

```python
from fastq_filter import filter_fastq

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

### Command-line Interface

`fastq_filter.py` can also be run directly from the command line:

```bash
# Basic usage
python fastq_filter.py --input seqs.fastq --output filtered_seqs.fastq

# With all filters and a custom log file
python fastq_filter.py --input seqs.fastq --output filtered_seqs.fastq \
    --length-bounds 25,150 --quality 20 --gc-bounds 40,60 \
    --log-file seqs_filter.log
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

**Note:** This parser is tailored to a specific output layout. For more robust parsing, consider using Biopython's `Bio.SearchIO` with tabular BLAST output (`-outfmt 6` or `7`).

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

## 4. Tests (`fastq_filter_tests.py`)

The test suite covers `filter_fastq` using `pytest`. Tests are grouped into four classes:

| Class                | What is tested                                              |
|----------------------|-------------------------------------------------------------|
| `TestLengthFilter`   | Upper bound, lower bound, single-int shorthand              |
| `TestQualityFilter`  | Quality threshold excludes low-quality reads                |
| `TestGCFilter`       | GC upper bound excludes GC-rich reads                       |
| `TestFileIO`         | Output file is created; empty output when nothing passes    |
| `TestErrorHandling`  | `FileNotFoundError` for missing input; `ValueError` for invalid bounds |

Run the tests with:

```bash
pytest fastq_filter_tests.py
```

Test output is logged to `tests.log` via `loguru`.

---

## Error Handling

- `filter_fastq` raises `ValueError` when bounds are logically invalid (e.g., min > max, or negative values).
- `filter_fastq` raises `FileNotFoundError` when the input file does not exist.
- Sequence classes validate the alphabet at instantiation and raise `ValueError` if invalid characters are found.
- `NucleicAcidSequence` cannot be instantiated directly and raises `NotImplementedError`.

---

## Dependencies

- **Python** == 3.13
- **Biopython** == 1.86 – required for `fastq_filter.py` and `bio_files_processor.py`
- **loguru** == 0.7.3 – required for `fastq_filter.py` and `fastq_filter_tests.py`
- **pytest** == 8.4.2 – required for running `fastq_filter_tests.py`

Install all at once:

```bash
pip install biopython loguru pytest
```

## License

This code is provided as-is for educational and research purposes. Feel free to adapt and reuse.
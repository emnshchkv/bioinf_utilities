# Bioinformatics Utilities

Personal collection of Python utilities for processing DNA and RNA sequences, filtering FASTQ files, and processing biological file formats.

## Installation

The script requires the following modules in the `scripts` directory:
- `dna_rna_tools`
- `filter_fastq_tools` 
- `assistant_utilities`

## Modules

### DNA/RNA Sequence Tools

The `run_dna_rna_tools()` function provides various sequence manipulation and analysis operations:

```python
from bioinf_utilities import run_dna_rna_tools

# Single sequence operation
result = run_dna_rna_tools("ATCG", "reverse")

# Multiple sequences operation  
results = run_dna_rna_tools("ATCG", "GCTA", "complement")
```

**Available Tools:**
- `is_rna` - Check if sequence is RNA
- `is_dna` - Check if sequence is DNA  
- `is_nucleic_acid` - Check if sequence is nucleic acid
- `transcribe` - Transcribe DNA to mRNA
- `reverse` - Reverse sequence
- `complement` - Get complementary sequence
- `reverse_complement` - Get reverse complementary sequence

**Parameters:**
- `*funnel`: Variable number of arguments where the last argument is the tool name and preceding arguments are sequences

**Returns:**
- Single result, type `str`, if one sequence provided
- List of results, type `list`, if multiple sequences provided

**Raises:**
- `ValueError` if not enough arguments provided
- `KeyError` if unknown tool specified

### FASTQ File Filtering

The `filter_fastq()` function filters FASTQ files based on GC-content, length, and quality scores, saving results to an output file:

```python
from bioinf_utilities import filter_fastq

# Filter with default bounds (reads from file, writes to 'output.fastq')
filter_fastq("input.fastq")

# Filter with custom parameters
filter_fastq(
    "input.fastq",
    gc_bounds=(35, 60),      # GC content between 35% and 60%
    length_bounds=(50, 150), # Length between 50 and 150 bp
    quality_threshold=20     # Minimum quality score of 20
)

# Single bound examples
filter_fastq("input.fastq", gc_bounds=60)        # Maximum 60% GC
filter_fastq("input.fastq", length_bounds=100)   # Maximum 100 bp length
```

**Parameters:**
- `fastq_path`: Absolute path to input FASTQ file
- `gc_bounds`: GC-content bounds as integer/float (max) or tuple (min, max). Default: (0, 100)
- `length_bounds`: Length bounds as integer (max) or tuple (min, max). Default: (0, 2^32)
- `quality_threshold`: Minimum quality score threshold. Default: 0

**Returns:**
- `None` - Results are written to 'output.fastq' file in the same directory as input

**Features:**
- Processes standard 4-line FASTQ format (identifier, sequence, separator, quality)
- Applies all filters sequentially: GC content → length → quality
- Maintains original FASTQ structure in output
- Handles large files through streaming processing

### Biological File Processing

The `bio_files_processor` module provides utilities for working with common biological file formats:

#### FASTA File Conversion

```python
from bio_files_processor import convert_multiline_fasta_to_oneline

# Convert multi-line FASTA to single-line format
convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
```

**Function:** `convert_multiline_fasta_to_oneline()`
- Converts FASTA files with multi-line sequences into single-line format
- **Parameters:** `input_fasta` (input file path), `output_fasta` (optional output path)
- **Returns:** `None`

#### BLAST Output Parsing

```python
from bio_files_processor import parse_blast_output

# Extract first protein from Description column for each BLAST query
parse_blast_output("blast_results.txt", "extracted_proteins.txt")
```

**Function:** `parse_blast_output()`
- Parses BLAST output files and extracts the first protein from the Description column for each query
- **Parameters:** `input_file` (BLAST output path), `output_file` (optional output path)
- **Returns:** `None`

#### Gene Neighborhood Extraction

```python
from bio_files_processor import select_genes_from_gbk_to_fasta

# Extract neighboring genes from GenBank file
select_genes_from_gbk_to_fasta(
    "genome.gbk",
    2,  # n_before
    2,  # n_after
    "gene1", "gene2",  # genes of interest
    output_fasta="neighbors.fasta"
)
```

**Function:** `select_genes_from_gbk_to_fasta()`
- Extracts amino acid sequences of genes adjacent to genes of interest from GenBank files
- Genes of interest themselves are not included in the output
- **Parameters:** 
  - `input_gbk`: Path to GenBank file
  - `n_before`, `n_after`: Number of genes before and after each gene of interest (default: 1)
  - `*genes`: Variable number of gene names to find neighbors for
  - `output_fasta`: Output FASTA file name (default: "output.fasta")
- **Returns:** `None`

## Error Handling

- Functions include comprehensive type hints for better development experience
- Input validation and error checking for common file format issues
- Clear error messages for malformed biological data
- Graceful handling of file I/O operations

## Dependencies

- Standard Python libraries: `os`, `re`
- No external dependencies required for core functionality

## File Format Support

- **FASTA**: Multi-line and single-line sequence formatting
- **FASTQ**: Quality score filtering and sequence validation
- **GenBank**: Gene annotation and protein sequence extraction
- **BLAST**: Tabular output parsing and result extraction

## Example Workflow

```python
from bioinf_utilities import filter_fastq
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    select_genes_from_gbk_to_fasta
)

# Quality filter FASTQ data
filter_fastq(
    "raw_sequences.fastq",
    gc_bounds=(40, 65),
    length_bounds=(75, 200),
    quality_threshold=25
)

# Convert FASTA format
convert_multiline_fasta_to_oneline("raw_sequences.fasta", "processed.fasta")

# Extract gene neighborhoods for functional analysis
select_genes_from_gbk_to_fasta(
    "ecoli.gbk",
    3, 2,
    "rpsA", "rplB",
    output_fasta="ribosomal_neighbors.fasta"
)
```

## Output Files

- **FASTQ filtering**: Creates `output.fastq` in same directory as input
- **FASTA conversion**: Creates specified output file or `{input_name}_output.fasta`
- **BLAST parsing**: Creates specified output file or `{input_name}_output_parse_blast.txt`
- **Gene extraction**: Creates specified FASTA file or `output.fasta`
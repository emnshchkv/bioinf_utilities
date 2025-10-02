# Bioinformatics Utilities

Personal collection of Python utilities for processing DNA and RNA sequences and filtering FASTQ files.

## Usage

The script requires the following modules in the `scripts` directory:
   - `is_rna`
   - `is_dna` 
   - `is_nucleic_acid`
   - `transcribe`
   - `reverse`
   - `complement`
   - `reverse_complement`
   - `gc_score`
   - `quality_score`
   - `interval`

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

### FASTQ Filtering

The `filter_fastq()` function filters FASTQ sequences based on GC content, length, and quality scores:

```python
from bioinf_utilities import filter_fastq

# Example FASTQ data structure
fastq_data = {
    "seq1": ("ATCGATCG", "IIIIIIII"),
    "seq2": ("GCTAGCTA", "JJJJJJJJ")
}

# Filter with default bounds
filtered = filter_fastq(fastq_data)

# Filter with custom parameters
filtered = filter_fastq(
    fastq_data,
    gc_bounds=(40, 60),
    length_bounds=(50, 150),
    quality_threshold=20
)
```

**Parameters:**
- `seqs`: Dictionary with sequence names as keys and tuples of (sequence, quality) as values
- `gc_bounds`: GC content bounds as integer (exact) or tuple (min, max). Default: (0, 100)
- `length_bounds`: Length bounds as integer (exact) or tuple (min, max). Default: (0, 2^32)
- `quality_threshold`: Minimum quality score threshold. Default: 0

**Returns:**
- Dictionary containing sequences in FASTQ-format that pass all filters

## Error Handling

- Invalid tool names will raise `KeyError` with available options
- Insufficient arguments will raise `ValueError`
- Functions include type hints for better development experience
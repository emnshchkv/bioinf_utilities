from  Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def _normalize_bounds(bounds: int|tuple) -> tuple:
    return bounds if isinstance(bounds, tuple) else (0, bounds)

def filter_fastq(
    input_path: str, 
    output_path: str, 
    length_bounds: int|tuple = (0, 2**32), 
    quality_threshold: int|float = 0, 
    gc_bounds: int|tuple = (0, 100)
    ) -> None:
    """
    Filters sequences from a FASTQ file based on length, average quality, and GC content, and writes the passing sequences to a new FASTQ file.

    Args:
        input_path (str): Path to the input FASTQ file.
        output_path (str): Path to save the filtered FASTQ file.
        length_bounds (int or tuple of int, optional): 
            Minimum and maximum sequence lengths. 
            If a single int is given, treated as maximum length, minimum is 0. Default is (0, 2**32).
        quality_threshold (int or float, optional): Minimum average Phred quality score required. Default is 0.
        gc_bounds (int or tuple of int, optional): 
            Minimum and maximum GC content (%) allowed. 
            If a single int is given, treated as maximum, minimum is 0. Default is (0, 100).

    Returns:
        None: Filtered sequences are written directly to output_path.
    """
    
    length_left, length_right = _normalize_bounds(length_bounds)
    gc_left, gc_right = _normalize_bounds(gc_bounds)
        
    with open(output_path, "w") as out_handle:
        
        for record in SeqIO.parse(input_path, "fastq"):
            length = len(record)
            qualities = record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities)
            gc = gc_fraction(record.seq) * 100

            if (
                length_left <= length <= length_right
                and avg_quality >= quality_threshold
                and gc_left <= gc <= gc_right
            ):
                SeqIO.write(record, out_handle, "fastq")

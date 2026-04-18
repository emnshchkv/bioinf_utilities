import argparse
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from loguru import logger


def setup_logging(log_file: str = "fastq_filter.log") -> None:
    """
    Configure loguru to write to both a file and the console.

    Args:
        log_file (str): Path to the log file. Default is 'fastq_filter.log'.
    """
    logger.remove()

    logger.add(
        sys.stderr,
        level="DEBUG",
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> [<level>{level}</level>] {message}",
        colorize=True,
    )

    logger.add(
        log_file,
        level="DEBUG",
        format="{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}",
        rotation="10 MB",
        retention="7 days",
        encoding="utf-8",
    )


def _normalize_bounds(bounds: int | tuple) -> tuple:
    """
    Normalize a bounds parameter to return a tuple of two values.

    If parameter is a tuple, it is returned as it is.
    If parameter is an integer, it is treated as the upper bound, and the lower bound is set to 0.

    Args:
        bounds (int or tuple of int): A single integer or a tuple representing (lower, upper) bounds.

    Returns:f
        tuple: A tuple of two integers (lower_bound, upper_bound).
    """
    return bounds if isinstance(bounds, tuple) else (0, bounds)


def filter_fastq(
    input_path: str,
    output_path: str,
    length_bounds: int | tuple = (0, 2**32),
    quality_threshold: int | float = 0,
    gc_bounds: int | tuple = (0, 100),
) -> None:
    """
    Filters sequences from a FASTQ file based on length, average quality, and GC content,
    and writes the passing sequences to a new FASTQ file.

    Args:
        input_path (str): Path to the input FASTQ file.
        output_path (str): Path to save the filtered FASTQ file.
        length_bounds (int or tuple of int, optional):
            Minimum and maximum sequence length.
            If a single int is given, treated as maximum length, minimum is 0. Default is (0, 2**32).
        quality_threshold (int or float, optional): Minimum average Phred quality score required. Default is 0.
        gc_bounds (int or tuple of int, optional):
            Minimum and maximum GC content (%) allowed.
            If a single int is given, treated as maximum, minimum is 0. Default is (0, 100).

    Returns:
        None: Filtered sequences are written directly to output_path.
    """

    length_left, length_right = _normalize_bounds(length_bounds)

    if length_left > length_right:
        raise ValueError("Min length cannot be greater than max length.")
    if length_left < 0 or length_right <= 0:
        raise ValueError("Only positive values allowed.")

    gc_left, gc_right = _normalize_bounds(gc_bounds)

    if gc_left > gc_right:
        raise ValueError("Min GC-content cannot be greater than max GC-content.")
    if gc_left < 0 or gc_right <= 0:
        raise ValueError("Only positive values allowed.")

    logger.info(
        "Start filtering: input='{}', output='{}', "
        "length=({}, {}), quality>={:.1f}, gc=({:.1f}, {:.1f})",
        input_path,
        output_path,
        length_left,
        length_right,
        quality_threshold,
        gc_left,
        gc_right,
    )

    total, passed = 0, 0

    with open(output_path, "w") as out_handle:

        for record in SeqIO.parse(input_path, "fastq"):
            total += 1
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
                passed += 1

    logger.info(
        "Filtering complete: {} / {} sequences passed ({:.1f}%)",
        passed,
        total,
        (passed / total * 100) if total else 0,
    )


def _parse_bounds(value: str) -> int | tuple:
    """
    Parse a bounds argument from the command line.
    Accepts either a single integer ("150") or a tuple ("25,150").

    Args:
        value (str): The raw string value from argparse.

    Returns:
        int or tuple of int: Parsed bounds.
    """
    value = value.strip()
    parts = [p.strip() for p in value.split(",")]

    if len(parts) == 1:
        return int(parts[0])
    elif len(parts) == 2:
        return (int(parts[0]), int(parts[1]))
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid bounds format: '{value}'. "
            "Expected a single integer or two comma separated integers, e.g. '150' or '25,150'."
        )


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="FASTQ-Filter - filters sequences based on length, average quality, and GC-content",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
                    Examples:
                    # Basic usage
                    python fastq_filter.py --input seqs.fastq --output filtered_seqs.fastq
                    
                    # Using extended settings
                    python fastq_filter.py --input seqs.fastq --output filtered_seqs.fastq\\
                        --length-bounds 25,150 --quality 20 --gc-bounds 40,50 \\
                        --log-file seqs_filter.log
                    """,
    )

    parser.add_argument(
        "--input", "-i", type=str, required=True, help="Input FASTQ file path"
    )

    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Output FASTQ file path"
    )

    parser.add_argument(
        "--length-bounds",
        type=_parse_bounds,
        default=(0, 2**32),
        help="Min and max sequence length (default: 0,2**32)",
    )
    parser.add_argument(
        "--quality",
        type=float,
        default=0,
        help="Min average Phred quality score required (default: 0)",
    )
    parser.add_argument(
        "--gc-bounds",
        type=_parse_bounds,
        default=(0, 100),
        help="Min and max GC content (%) allowed (default: 0,100)",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default="fastq_filter.log",
        help="Path to the log file (default: fastq_filter.log)",
    )

    args = parser.parse_args()
    setup_logging(args.log_file)

    try:
        filter_fastq(
            input_path=args.input,
            output_path=args.output,
            length_bounds=args.length_bounds,
            quality_threshold=args.quality,
            gc_bounds=args.gc_bounds,
        )
    except FileNotFoundError as e:
        logger.error("Input file not found: {}", e)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error: {}", e)
        sys.exit(1)


if __name__ == "__main__":
    main()

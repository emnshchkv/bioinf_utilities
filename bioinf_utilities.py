import os
from scripts import filter_fastq_tools, assistant_ulitities
from typing import Union

def filter_fastq(
    fastq_path: str,
    gc_bounds: int | float | tuple[int | float, int | float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int | float = 0,
) -> None:
    """
    Filters *.fastq with manual settings and save it to 'output.fastq'.

    Agruments:
    fastq_path: absolute path to *.fastq
    gc_bounds: int | tuple
    length_bounds: int | tuple
    quality_threshold: int

    Returns:
    None
    """
    (gc_lower_bound, gc_upper_bound), (len_lower_bound, len_upper_bound) = (
        assistant_ulitities.make_interval(gc_bounds),
        assistant_ulitities.make_interval(length_bounds),
    )
    input_path = os.path.join(fastq_path)
    output_path = assistant_ulitities.make_output_path(fastq_path, "output.fastq")

    with open(input_path, "r") as input_fastq, open(output_path, "a") as output_fastq:
        while True:
            identifier = input_fastq.readline().strip()
            if not identifier:
                break
            sequence = input_fastq.readline().strip()
            _ = input_fastq.readline().strip()
            quality = input_fastq.readline().strip()
            if filter_fastq_tools.is_suitable_fastq(
                sequence,
                quality,
                gc_lower_bound,
                gc_upper_bound,
                len_lower_bound,
                len_upper_bound,
                quality_threshold,
            ):
                output_fastq.write(
                    identifier + "\n" + sequence + "\n" + _ + "\n" + quality
                )

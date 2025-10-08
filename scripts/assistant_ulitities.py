from typing import Union
import os


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

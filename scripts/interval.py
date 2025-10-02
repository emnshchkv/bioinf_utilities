from typing import Union


def interval(
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

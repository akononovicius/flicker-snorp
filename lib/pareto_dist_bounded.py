import numpy as np


def sample(
    power: float,
    low: float = 1,
    high: float = 1000,
    size: int | tuple = 1,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Sample from the (un)bounded Pareto distribution.

    Input:
        power:
            Power law exponent parameter. Note that the actual value
            of the exponent in the PDF is `power+1`.
        low: (default: 1)
            Location of the lower truncation.
        high: (default: 1000)
            Location of the upper truncation.
        size: (default: 1)
            Size of the desired sample.
        rng: (default: None)
            Numpy RNG object. If not specified, then a default RNG
            instance will be created within this function.

    Output:
        Samples arranged into Numpy array of predetermined size, or
        singular real value.
    """
    if rng is None:
        rng = np.random.default_rng()

    if high == np.inf:
        if size == 1:
            return np.array((rng.pareto(power) + 1) * low)
        else:
            return (rng.pareto(power, size=size) + 1) * low

    if size == 1:
        scale = high / low
        u = rng.random()
        result = (1 - u + u / (scale**power)) ** (-1 / power)
        result = result * low
    else:
        n_elems = int(np.prod(size))
        result = np.array(
            [sample(power, low=low, high=high, size=1, rng=rng) for _ in range(n_elems)]
        )
        result.shape = size

    return result

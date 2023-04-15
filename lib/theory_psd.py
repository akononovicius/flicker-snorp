import numpy as np
from mpmath import gamma as mpmath_gamma  # type: ignore

euler_gamma = 0.5772156649015328606065120


def __gamma(x: float) -> float:
    return float(mpmath_gamma(x))


def get_poiss_poiss_psd(
    freqs: np.ndarray,
    pulse_magnitude: float,
    mean_pulse: float,
    mean_gap: float,
) -> np.ndarray:
    """Calculate theoretical PSD for Poissonian pulses and gaps.

    Input:
        freqs:
            Desired frequencies.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        mean_pulse:
            Average pulse duration.
        mean_gap:
            Average gap duration.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    angular_freqs = 2 * np.pi * freqs

    gamma_theta = 1 / mean_pulse
    gamma_tau = 1 / mean_gap

    nu_bar = 1 / (mean_pulse + mean_gap)

    const_term = 4 * (pulse_magnitude**2) * nu_bar

    return const_term / ((gamma_theta + gamma_tau) ** 2 + angular_freqs**2)


def get_long_poiss_bounded_pareto_psd(
    freqs: np.ndarray,
    pulse_magnitude: float,
    mean_pulse: float,
    min_gap: float,
    max_gap: float,
    power_gap: float,
) -> np.ndarray:
    """Calculate theoretical PSD for long Poissonian pulses and bounded Pareto gaps.

    Input:
        freqs:
            Desired frequencies.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        mean_pulse:
            Average pulse duration.
        min_gap:
            Minimum gap duration.
        max_gap:
            Maximum gap duration.
        power_gap:
            Power law exponent of the bounded Pareto gap distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    angular_freqs = 2 * np.pi * freqs

    mean_gap = __get_pareto_mean(min_gap, max_gap, power_gap)
    nu_bar = 1 / (mean_pulse + mean_gap)

    result = np.array([])

    if power_gap == 1:
        result = min_gap / freqs
    elif power_gap < 2:
        const_term = 4 * (min_gap**2)
        const_term = const_term * __gamma(1 - power_gap)
        const_term = const_term * np.cos(np.pi * power_gap / 2)
        result = const_term * ((angular_freqs * min_gap) ** (power_gap - 2))
    else:
        gap_ratio = power_gap / (power_gap - 2)
        result = np.zeros(freqs.shape) + 2 * (min_gap**2) * gap_ratio

    return (pulse_magnitude**2) * nu_bar * result


def get_short_poiss_bounded_pareto_psd(
    freqs: np.ndarray,
    pulse_magnitude: float,
    mean_pulse: float,
    min_gap: float,
    max_gap: float,
    power_gap: float,
) -> np.ndarray:
    """Calculate theoretical PSD for short Poissonian pulses and bounded Pareto gaps.

    Input:
        freqs:
            Desired frequencies.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        mean_pulse:
            Average pulse duration.
        min_gap:
            Minimum gap duration.
        max_gap:
            Maximum gap duration.
        power_gap:
            Power law exponent of the bounded Pareto gap distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    angular_freqs = 2 * np.pi * freqs

    mean_gap = __get_pareto_mean(min_gap, max_gap, power_gap)
    nu_bar = 1 / (mean_pulse + mean_gap)

    result = np.array([])
    const_term = 0.0

    if power_gap == 1:
        result = 1 / freqs
        result = result / min_gap
        result = result / (
            (np.pi**2)
            + 4 * ((1 - euler_gamma - np.log(angular_freqs * min_gap)) ** 2)
        )
    elif power_gap < 1:
        const_term = np.cos(np.pi * power_gap / 2) / __gamma(1 - power_gap)
        result = const_term / ((angular_freqs * min_gap) ** power_gap)
    elif power_gap < 2:
        const_term = ((power_gap - 1) / power_gap) ** 2
        const_term = const_term * np.cos(np.pi * power_gap / 2) * __gamma(1 - power_gap)
        result = const_term / ((angular_freqs * min_gap) ** (2 - power_gap))
    else:
        const_term = (power_gap - 1) ** 2
        const_term = const_term / (2 * (power_gap - 2) * power_gap)
        result = np.zeros(freqs.shape) + const_term

    return 4 * (pulse_magnitude**2) * nu_bar * (mean_pulse**2) * result


def __get_pareto_mean(
    low: float,
    high: float,
    power: float,
) -> float:
    """Calculcate average value of the bounded Pareto distribution.

    Input:
        low:
            Lower truncation point.
        high:
            Upper truncation point.
        power:
            Power law exponent of the distribution. Note that the
            actual value of the exponent in the PDF is `power+1`.

    Output:
        Average value.
    """
    if power == 1:
        return high * low / (high - low) * np.log(high / low)
    term_1 = (low**power) / (1 - (low / high) ** power)
    term_2 = power / (power - 1)
    term_3 = 1 / (low ** (power - 1)) - 1 / (high ** (power - 1))
    return term_1 * term_2 * term_3


def get_any_bounded_pareto_psd(
    freqs: np.ndarray,
    nu_bar: float,
    pulse_magnitude: float,
    min_gap: float,
    max_gap: float,
    power_gap: float,
):
    """Calculate theoretical PSD for any pulses and bounded Pareto gaps.

    Input:
        freqs:
            Desired frequencies.
        nu_bar:
            Average number of pulses per unit time.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        min_gap:
            Minimum gap duration.
        max_gap:
            Maximum gap duration.
        power_gap:
            Power law exponent of the bounded Pareto gap distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    if power_gap != 1:
        result = np.empty(freqs.shape)
        result[:] = np.nan
        return result

    return (pulse_magnitude**2) * nu_bar * min_gap / freqs


def get_const_bounded_pareto_psd(
    freqs: np.ndarray,
    pulse_magnitude: float,
    fixed_pulse: float,
    min_gap: float,
    max_gap: float,
    power_gap: float,
) -> np.ndarray:
    """Calculate theoretical PSD for fixed pulses and bounded Pareto gaps.

    Input:
        freqs:
            Desired frequencies.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        fixed_pulse:
            Fixed pulse duration.
        min_gap:
            Minimum gap duration.
        max_gap:
            Maximum gap duration.
        power_gap:
            Power law exponent of the bounded Pareto gap distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    mean_gap = __get_pareto_mean(min_gap, max_gap, power_gap)
    nu_bar = 1 / (fixed_pulse + mean_gap)

    return get_any_bounded_pareto_psd(
        freqs, nu_bar, pulse_magnitude, min_gap, max_gap, power_gap
    )


def get_uniform_bounded_pareto_psd(
    freqs: np.ndarray,
    pulse_magnitude: float,
    min_pulse: float,
    max_pulse: float,
    min_gap: float,
    max_gap: float,
    power_gap: float,
) -> np.ndarray:
    """Calculate theoretical PSD for uniform pulses and bounded Pareto gaps.

    Input:
        freqs:
            Desired frequencies.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        min_pulse:
            Minimum pulse duration.
        max_pulse:
            Maximum pulse duration.
        min_gap:
            Minimum gap duration.
        max_gap:
            Maximum gap duration.
        power_gap:
            Power law exponent of the bounded Pareto gap distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    mean_gap = __get_pareto_mean(min_gap, max_gap, power_gap)
    mean_uniform = (min_pulse + max_pulse) / 2
    nu_bar = 1 / (mean_uniform + mean_gap)

    return get_any_bounded_pareto_psd(
        freqs, nu_bar, pulse_magnitude, min_gap, max_gap, power_gap
    )


def get_double_bounded_pareto_psd(
    freqs: np.ndarray,
    pulse_magnitude: float,
    min_pulse: float,
    max_pulse: float,
    power_pulse: float,
    min_gap: float,
    max_gap: float,
    power_gap: float,
) -> np.ndarray:
    """Calculate theoretical PSD for bounded Pareto pulses and gaps.

    Input:
        freqs:
            Desired frequencies.
        pulse_magnitude:
            Fixed magnitude of the pulses in the signal.
        min_pulse:
            Minimum pulse duration.
        max_pulse:
            Maximum pulse duration.
        power_pulse:
            Power law exponent of the bounded Pareto pulse distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.
        min_gap:
            Minimum gap duration.
        max_gap:
            Maximum gap duration.
        power_gap:
            Power law exponent of the bounded Pareto gap distribution. Note
            that the actual value of the exponent in the PDF is `power+1`.

    Output:
        A theoretical estimate of PSD values at the desired
        frequencies.
    """
    mean_gap = __get_pareto_mean(min_gap, max_gap, power_gap)
    mean_pulse = __get_pareto_mean(min_pulse, max_pulse, power_pulse)
    nu_bar = 1 / (mean_pulse + mean_gap)

    return get_any_bounded_pareto_psd(
        freqs, nu_bar, pulse_magnitude, min_gap, max_gap, power_gap
    )

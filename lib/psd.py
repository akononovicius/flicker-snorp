import numpy as np


def get_snorp_psd(
    freqs: np.ndarray,
    pulse_durations: np.ndarray,
    gap_durations: np.ndarray,
    pulse_magnitude: float = 1,
) -> np.ndarray:
    """Calculate PSD of the signal with non-overlapping rectangular pulses.

    Input:
        freqs:
            Frequencies for which to obtain the estimates
            of the PSD.
        pulse_durations:
            List containing durations of each pulse in
            the signal.
        gap_durations:
            List containing durations of each interpulse
            gap in the signal
        pulse_magnitude: (default: 1)
            Fixed magnitude of the pulses in the signal.

    Output:
        An estimate of the PSD at given frequencies.

    Note: This function is applicable only when the pulses
        are of rectangular shape and have the same magnitude.
        By default unit magnitude is assumed.
    """
    angular_freqs = 2 * np.pi * freqs

    # our simplification of the Fourier transform formula requires having not
    # only pulse or gap durations, but also the time moment when the respective
    # pulses or gaps have started
    pulse_starts = (
        np.cumsum(pulse_durations) + np.cumsum(gap_durations) - pulse_durations
    )
    gap_starts = pulse_starts - gap_durations

    # to avoid artificats for the lowest frequences we need to subtact the mean
    # magnitude of the signal from the series prior to applying Fourier
    # transform
    total_duration = pulse_starts[-1] + pulse_durations[-1]
    in_pulse_time = np.sum(pulse_durations)
    mean_magnitude = pulse_magnitude * in_pulse_time / total_duration
    adjusted_pulse_magnitude = pulse_magnitude - mean_magnitude
    adjusted_gap_magnitude = -mean_magnitude

    normalization = 2 / total_duration

    return normalization * np.array(
        [
            __get_snorp_psd(
                omega,
                pulse_durations,
                pulse_starts,
                adjusted_pulse_magnitude,
                gap_durations,
                gap_starts,
                adjusted_gap_magnitude,
            )
            for omega in angular_freqs
        ]
    )


def __get_snorp_psd(
    angular_freq: float,
    pulse_durations: np.ndarray,
    pulse_starts: np.ndarray,
    pulse_magnitude: float,
    gap_durations: np.ndarray,
    gap_starts: np.ndarray,
    gap_magnitude: float,
) -> float:
    """Calculate non-normalized PSD of the signal non-overlapping rectangular pulses.

    Input:
        angular_freq:
            Angular frequency for which to obtain the estimate
            of the PSD.
        pulse_durations:
            List containing durations of each pulse in the
            signal.
        pulse_starts:
            List containing start times of each pulse in the
            signal.
        pulse_magnitude:
            Effective magnitude of the pulses in the signal.
        gap_durations:
            List containing durations of each interpulse gap
            in the signal.
        gap_starts:
            List containing start times of each interpulse
            gap in the signal.
        gap_magnitude:
            Effective magnitude of the interpulse gaps in
            the signal

    Output:
        An estimate of the non-normalized PSD at a given angular
        frequency.
    """
    fourier = __get_rect_fourier(
        angular_freq, pulse_magnitude, pulse_durations, pulse_starts
    ) + __get_rect_fourier(angular_freq, gap_magnitude, gap_durations, gap_starts)
    return np.real(fourier) ** 2 + np.imag(fourier) ** 2


def __get_rect_fourier(
    angular_freq: float,
    magnitude: float,
    durations: np.ndarray,
    starts: np.ndarray,
) -> float:
    """Calculate Fourier transform of a rectangular pulse.

    Input:
        angular_freq:
            Desired angular frequency.
        magnitude:
            Magnitude of the pulses.
        durations:
            Durations of the pulses.
        starts:
            Start times of the pulses.

    Output:
        A Fourier transform (complex number) at a given angular
        frequency.

    Note: This function is applicable only when the pulses
        are of rectangular. This is encoded in the analytically
        precomputed expression for the "profile".
    """
    constant_terms = magnitude * (1j / angular_freq)
    profiles = np.exp(-1j * angular_freq * durations) - 1
    variable_terms = np.exp(-1j * angular_freq * starts) * profiles
    return constant_terms * np.sum(variable_terms)

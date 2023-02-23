import numpy as np
import typer

from lib.pareto_dist_bounded import sample
from lib.psd import (get_long_poiss_bounded_pareto_psd,
                     get_short_poiss_bounded_pareto_psd, get_snorp_psd)


def simulate_duration(
    T,
    pulse_magnitude,
    mean_pulse,
    min_gap,
    max_gap,
    power_gap,
    rng,
):
    # sample distributions
    sample_pulse = rng.exponential

    def sample_gap(power, low=1, high=1000, size=1):
        return sample(power, low=low, high=high, size=size, rng=rng)

    # main
    t_pulse = 0
    t_gap = 0
    pulse_duration = []
    gap_duration = []
    while t_pulse + t_gap < T:
        gap = sample_gap(power_gap, low=min_gap, high=max_gap)
        pulse = sample_pulse(mean_pulse)
        if t_gap + t_pulse + gap > T:
            gap = T - t_pulse - t_gap
            pulse = 0
        t_gap = t_gap + gap
        if pulse > 0 and t_gap + t_pulse + pulse > T:
            pulse = T - t_pulse - t_gap
        t_pulse = t_pulse + pulse
        pulse_duration += [pulse]
        gap_duration += [gap]
    return np.array(pulse_duration), np.array(gap_duration)


def main(
    repeats: int = 1,
    duration: float = 1e6,
    pulse_magnitude: float = 1,
    mean_pulse: float = 1,
    min_gap: float = 1,
    max_gap: float = 1e4,
    power_gap: float = 0.5,
    min_freq: float = -1,
    max_freq: float = -1,
    n_freq: int = 100,
    archive_dir: str = "data",
    seed: int = -1,
):
    """Simulate rectangular SNORP with Poissonian pulses and bounded Pareto gaps.

    Input:
        repeats: (default: 1)
            Number of times to generate SNORP. Resulting PSD
            will be averaged over all runs.
        duration: (default: 1e6)
            Duration over which to simulate.
        pulse_magnitude: (default: 1)
            Fixed magnitude of the pulses in the signal.
        mean_pulse: (default: 1)
            Mean pulse duration.
        min_gap: (default: 1)
            Minimum gap duration.
        max_gap: (default: 1e4)
            Maximum gap duration. If negative number is
            passed, then infinite value is used instead.
        power_gap: (default: 0.5)
            Power law exponent of gap duration distribution.
        min_freq: (default: -1)
            Minimum frequency to observe. If negative value is
            passed (which is the default), then the minimum
            frequency is calculated automatically based on
            simulation parameters.
        max_freq: (default: -1)
            Maximum frequency to observe. If negative value is
            passed (which is the default), then the maximum
            frequency is calculated automatically based on
            simulation parameters.
        n_freq: (default: 100)
            Number of frequencies to consider within the given
            (or automatically selected) range. Includes end
            points.
        archive_dir (default: "data")
            Folder in which to save output files.
        seed: (default: -1)
            RNG seed. If negative value is passed (which is the
            default), then seed will be randomly generated by
            `np.random.randint(0, int(2**20))`

    Output:
        Function returns nothing, but saves one file, which
        contains the numerically calculated PSD and its
        theoretical estimate.
    """
    # auto-generate seed
    if seed < 0:
        np.random.seed()
        seed = np.random.randint(0, int(2**20))

    # RNG setup
    rng = np.random.default_rng(seed)

    # simulation archival setup
    model_info = f"poiss{mean_pulse*10000:.0f}.pareto{power_gap*100:.0f}_{min_gap*1000:.0f}_{max_gap:.0f}"
    simulation_filename = f"{model_info}.seed{seed:d}"
    psd_path = f"{archive_dir}/{simulation_filename}.psd.csv"

    # negative max_gap implies inf max_gap
    max_gap_eff = max_gap
    if max_gap < 0:
        max_gap = np.inf
        max_gap_eff = duration

    # set frequency range
    if max_freq < 0:
        max_freq = (1 / max_gap_eff) * 0.1 / (2 * np.pi)
    if min_freq < 0:
        min_freq = (1 / min_gap) * 10 / (2 * np.pi)
    freqs = np.logspace(np.log10(min_freq), np.log10(max_freq), n_freq)
    freqs = np.unique(np.round(duration * freqs)) / duration  # only natural freqs
    n_freq = len(freqs)

    # main simulation loop
    sim_psds = np.zeros((repeats, n_freq))
    for sim_idx in range(repeats):
        pulse_duration, gap_duration = simulate_duration(
            duration,
            pulse_magnitude,
            mean_pulse,
            min_gap,
            max_gap,
            power_gap,
            rng,
        )
        sim_psds[sim_idx, :] = get_snorp_psd(
            freqs,
            pulse_duration,
            gap_duration,
            pulse_magnitude=pulse_magnitude,
        )

    # numerical PSD
    sim_psd = np.mean(sim_psds, axis=0)

    # theoretical PSD
    if mean_pulse < min_gap:
        theory_psd = get_short_poiss_bounded_pareto_psd(
            freqs,
            pulse_magnitude,
            mean_pulse,
            min_gap,
            max_gap_eff,
            power_gap,
        )
    else:
        theory_psd = get_long_poiss_bounded_pareto_psd(
            freqs,
            pulse_magnitude,
            mean_pulse,
            min_gap,
            max_gap_eff,
            power_gap,
        )

    np.savetxt(
        psd_path,
        np.log10(np.vstack((freqs, sim_psd, theory_psd)).T),
        delimiter=",",
        fmt="%.4f",
    )


if __name__ == "__main__":
    typer.run(main)

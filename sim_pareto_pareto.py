import numpy as np
import typer

from lib.pareto_dist_bounded import sample
from lib.psd import get_double_bounded_pareto_psd, get_snorp_psd


def main(
    repeats: int = 1,
    n_events: int = 10000,
    pulse_magnitude: float = 1,
    min_pulse: float = 1,
    max_pulse: float = 1e4,
    power_pulse: float = 1,
    min_gap: float = 1,
    max_gap: float = 1e4,
    power_gap: float = 1,
    min_freq: float = -1,
    max_freq: float = -1,
    n_freq: int = 100,
    archive_dir: str = "data",
    seed: int = -1,
):
    """Simulate rectangular SNORP and bounded Pareto pulses and gaps.

    Input:
        repeats: (default: 1)
            Number of times to generate SNORP. Resulting PSD
            will be averaged over all runs.
        n_events: (default: 10000)
            Number of pulses to simulate.
        pulse_magnitude: (default: 1)
            Fixed magnitude of the pulses in the signal.
        min_pulse: (default: 1)
            Minimum pulse duration.
        max_pulse: (default: 1e4)
            Maximum pulse duration.
        power_pulse: (default: 1)
            Power law exponent of pulse duration distribution.
        min_gap: (default: 1)
            Minimum gap duration.
        max_gap: (default: 1e4)
            Maximum gap duration.
        power_gap: (default: 1)
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
    model_info = f"pareto{power_pulse*100:.0f}_{min_pulse*1000:.0f}_{max_pulse:.0f}.pareto{power_gap*100:.0f}_{min_gap*1000:.0f}_{max_gap:.0f}"
    simulation_filename = f"{model_info}.seed{seed:d}"
    psd_path = f"{archive_dir}/{simulation_filename}.psd.csv"

    # set frequency range
    if max_freq < 0:
        max_freq = (1 / np.max([max_gap, max_pulse])) * 0.1 / (2 * np.pi)
    if min_freq < 0:
        min_freq = (1 / np.min([min_gap, min_pulse])) * 10 / (2 * np.pi)
    freqs = np.logspace(np.log10(min_freq), np.log10(max_freq), n_freq)

    # main simulation loop
    sim_psds = np.zeros((repeats, n_freq))
    for sim_idx in range(repeats):
        pulse_duration = sample(
            power_pulse, low=min_pulse, high=max_pulse, size=n_events, rng=rng
        )
        gap_duration = sample(
            power_gap, low=min_gap, high=max_gap, size=n_events, rng=rng
        )
        sim_psds[sim_idx, :] = get_snorp_psd(
            freqs,
            pulse_duration,  # type: ignore
            gap_duration,  # type: ignore
            pulse_magnitude=pulse_magnitude,
        )

    # numerical PSD
    sim_psd = np.mean(sim_psds, axis=0)

    # theoretical PSD
    theory_psd = get_double_bounded_pareto_psd(
        freqs,
        pulse_magnitude,
        min_pulse,
        max_pulse,
        power_pulse,
        min_gap,
        max_gap,
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

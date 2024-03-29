import numpy as np
import typer

from lib.psd import get_snorp_psd
from lib.theory_psd import get_poiss_poiss_psd


def main(
    repeats: int = 1,
    n_events: int = 10000,
    pulse_magnitude: float = 1,
    mean_pulse: float = 1,
    mean_gap: float = 1,
    min_freq: float = -1,
    max_freq: float = -1,
    n_freq: int = 100,
    archive_dir: str = "data",
    seed: int = -1,
) -> None:
    """Simulate rectangular SNORP with Poissonian durations.

    Input:
        repeats: (default: 1)
            Number of times to generate SNORP. Resulting PSD
            will be averaged over all runs.
        n_events: (default: 10000)
            Number of pulses to simulate.
        pulse_magnitude: (default: 1)
            Fixed magnitude of the pulses in the signal.
        mean_pulse: (default: 1)
            Mean pulse duration.
        mean_gap: (default: 1)
            Mean gap duration.
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

    # sample distributions
    sample_pulse = rng.exponential
    sample_gap = rng.exponential

    # simulation archival setup
    model_info = f"poiss{mean_pulse*10:.0f}.poiss{mean_gap*10:.0f}"
    simulation_filename = f"{model_info}.seed{seed:d}"
    psd_path = f"{archive_dir}/{simulation_filename}.psd.csv"

    # set frequency range
    if min_freq < 0:
        min_freq = 1 / (n_events * (mean_gap + mean_pulse))
    if max_freq < 0:
        max_freq = 2 * n_events / np.min([mean_gap, mean_pulse])
    freqs = np.logspace(np.log10(min_freq), np.log10(max_freq), n_freq)

    # main simulation loop
    sim_psds = np.zeros((repeats, n_freq))
    for sim_idx in range(repeats):
        pulse_duration = sample_pulse(scale=mean_pulse, size=n_events)
        gap_duration = sample_gap(scale=mean_gap, size=n_events)
        sim_psds[sim_idx, :] = get_snorp_psd(
            freqs,
            pulse_duration,
            gap_duration,
            pulse_magnitude=pulse_magnitude,
        )

    # numerical PSD
    sim_psd = np.mean(sim_psds, axis=0)

    # theoretical PSD
    theory_psd = get_poiss_poiss_psd(freqs, pulse_magnitude, mean_pulse, mean_gap)

    np.savetxt(
        psd_path,
        np.log10(np.vstack((freqs, sim_psd, theory_psd)).T),
        delimiter=",",
        fmt="%.4f",
    )


if __name__ == "__main__":
    typer.run(main)

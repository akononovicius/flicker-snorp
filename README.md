# 1/f (flicker) noise from the sequence of non-overlapping rectangular pulses

This repository replicates numerical simulation conducted in [1]. Figure
generation scripts are also included.

Main runable files are "sim.sh", which conducts simulation, and "fig.sh",
which produces figures seen in [1]. Model specific simulation files are
named `sim_<model>.py` (where `<model>` is usually an abbreviation of
distributions from which pulse and gap duration are sampled). Figure
specific files are named `fig_<figure>.py` (where `<figure>` matches
corresponding figure name).

## References

1. A. Kononovicius, B. Kaulakys. *1/f noise from the sequence of
   nonoverlapping rectangular pulses.* Physical Review E **107**: 034117
   (2023). [doi:
   10.1103/PhysRevE.107.034117](https://doi.org/10.1103/PhysRevE.107.034117).
   [arXiv:2210.11792
   [cond-mat.stat-mech]](https://arxiv.org/abs/2210.11792).

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

plt.rc("text", usetex=True)


fNs = np.array(
    [
        "data/poiss1000000.pareto100_1000_-1.seed13007.psd.csv",
        "data/poiss1000000.pareto100_1000_-1.seed22675.psd.csv",
        "data/poiss1000000.pareto100_1000_-1.seed1081.psd.csv",
        "data/poiss1000000.pareto100_1000_-1.seed4859.psd.csv",
    ]
)


with PdfPages("figs/nonergodic.pdf") as pdfFile:
    fig = plt.figure(figsize=(3, 2))
    plt.loglog()
    plt.ylim([3e-6, 3e8])
    plt.xlim([3e-11, 3e0])

    ax = plt.gca()
    ax.grid(visible=True, linestyle="--", dashes=(15, 15), c="#cccccc", lw=0.3)
    ax.get_yaxis().set_tick_params(direction="out")
    ax.get_xaxis().set_tick_params(direction="out")
    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$S(f)$")

    ax.set_xticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    ax.set_yticks([1e-5, 1e-2, 1e1, 1e4, 1e7])

    data = 10 ** np.loadtxt(fNs[0], delimiter=",")
    plt.plot(data[::4, 0], data[::4, 1], "md", lw=4)
    data = 10 ** np.loadtxt(fNs[1], delimiter=",")
    plt.plot(data[::8, 0], data[::8, 1], "bs", lw=4)
    data = 10 ** np.loadtxt(fNs[2], delimiter=",")
    plt.plot(data[::10, 0], data[::10, 1], "g^", lw=4)
    data = 10 ** np.loadtxt(fNs[3], delimiter=",")
    plt.plot(data[::14, 0], data[::14, 1], "ro", lw=4)
    data = 10 ** np.loadtxt(fNs[0], delimiter=",")
    plt.plot(data[:, 0], 8e-3 / data[:, 0], "k--")

    plt.tight_layout(w_pad=1, h_pad=1, pad=0.5)

    pdfFile.savefig(fig)
    plt.close()

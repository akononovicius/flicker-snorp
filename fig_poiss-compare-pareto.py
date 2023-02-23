import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

plt.rc("text", usetex=True)


fNs = np.array(
    [
        "data/poiss10.pareto100_1000_10000.seed3981.psd.csv",
        "data/poiss10000000.pareto100_1000_10000.seed1186.psd.csv",
    ]
)


with PdfPages("figs/poiss-compare-pareto.pdf") as pdfFile:
    fig = plt.figure(figsize=(3, 2))
    plt.loglog()
    plt.ylim([3e-5, 3e2])
    plt.xlim([3e-6, 3e0])

    ax = plt.gca()
    ax.grid(False)
    # ax.grid(visible=True, linestyle="--", dashes=(15, 15), c="#cccccc", lw=0.3)
    ax.get_yaxis().set_tick_params(direction="out")
    ax.get_xaxis().set_tick_params(direction="out")
    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$S(f)$")

    ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0])
    ax.set_yticks([1e-4, 1e-2, 1e0, 1e2])

    data = 10 ** np.loadtxt(fNs[0], delimiter=",")
    data[:, 1] = data[:, 1]
    plt.plot(data[:, 0], data[:, 1], "g-", lw=4)
    plt.plot(data[:, 0], data[:, 2], "k--")
    data = 10 ** np.loadtxt(fNs[1], delimiter=",")
    plt.plot(data[:, 0], data[:, 1], "r-", lw=4)
    plt.plot(data[:, 0], data[:, 2], "k-")

    plt.tight_layout(w_pad=1, h_pad=1, pad=0.5)

    pdfFile.savefig(fig)
    plt.close()

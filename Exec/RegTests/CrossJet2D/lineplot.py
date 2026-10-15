"""
Plot over line along the lower wall (y = 0, first cell row of the finest
level) around the injector, for the 2D CrossJet runs.

    python3 lineplot.py RUN_DIR [RUN_DIR ...] [--times 10 30 50 70] [--out name.png]

With one run directory: one curve per time. With several: the runs are
overlaid at each time (one colour per time, one line style per run), so
pure PeleC and CLEM can be compared directly.
"""
import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

yt.set_log_level(50)

X_LO, X_HI = 3.9, 5.8  # cm, the refined box around the injector
JET_X, JET_D = 4.2, 0.1

FIELDS = [
    ("pressure", "p [dyn/cm$^2$]"),
    ("Temp", "T [K]"),
    ("density", r"$\rho$ [g/cm$^3$]"),
    ("Y(H2)", r"$Y_{H_2}$"),
    ("MachNumber", "Mach"),
    ("x_velocity", "u [cm/s]"),
]


def plotfile_at(run, t_us):
    """The plotfile of run closest to time t_us (microseconds)."""
    best, best_dt = None, np.inf
    for f in sorted(glob.glob(os.path.join(run, "plt[0-9]*"))):
        with open(os.path.join(f, "Header")) as h:
            lines = h.read().split("\n")
        n_var = int(lines[1])
        t = float(lines[2 + n_var + 1]) * 1.0e6
        if abs(t - t_us) < best_dt:
            best, best_dt = (f, t), abs(t - t_us)
    return best


def wall_line(pltfile):
    """x and the fields on the first row of cells above y = 0, finest level."""
    ds = yt.load(pltfile)
    lev = ds.index.max_level
    dx = float(ds.index.get_smallest_dx())
    i_lo = int(np.floor(X_LO / dx))
    i_hi = int(np.ceil(X_HI / dx))
    cg = ds.covering_grid(lev, left_edge=[i_lo * dx, 0.0, 0.0],
                          dims=[i_hi - i_lo, 1, 1])
    x = (np.arange(i_lo, i_hi) + 0.5) * dx
    data = {name: np.asarray(cg["boxlib", name])[:, 0, 0] for name, _ in FIELDS}
    return x, data, dx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--times", nargs="+", type=float, default=[10, 30, 50, 70])
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    styles = ["-", "--", ":", "-."]
    colours = plt.cm.viridis(np.linspace(0.0, 0.85, len(args.times)))
    fig, axes = plt.subplots(len(FIELDS), 1, figsize=(9, 2.3 * len(FIELDS)), sharex=True)

    for ir, run in enumerate(args.runs):
        label_run = os.path.basename(os.path.normpath(run))
        for it, t_us in enumerate(args.times):
            found = plotfile_at(run, t_us)
            if found is None:
                continue
            pltfile, t = found
            x, data, dx = wall_line(pltfile)
            for ax, (name, _) in zip(axes, FIELDS):
                label = f"t = {t:.0f} $\\mu$s" if len(args.runs) == 1 else f"{label_run}, {t:.0f} $\\mu$s"
                ax.plot(x, data[name], styles[ir % len(styles)], color=colours[it], lw=1.2,
                        label=label)

    for ax, (_, ylabel) in zip(axes, FIELDS):
        ax.axvspan(JET_X - 0.5 * JET_D, JET_X + 0.5 * JET_D, color="0.85", zorder=0)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    axes[0].set_title(f"Lower wall, first cell row (y = {0.5 * dx * 10:.3f} mm), injector slot shaded")
    axes[0].legend(fontsize=7, ncol=2)
    axes[-1].set_xlabel("x [cm]")
    fig.tight_layout()

    out = args.out or os.path.join(args.runs[0], "wall_line.png")
    fig.savefig(out, dpi=150)
    print(out)


if __name__ == "__main__":
    main()

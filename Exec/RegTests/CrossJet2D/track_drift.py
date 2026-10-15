"""
Energy drift between the CLEM lines and the LES over time, from run.log.

    python3 track_drift.py RUN_DIR [--out drift.png]

Per level and per step, the CLEM step prints (advance_clem.cpp):
  - "line energy anchored to the LES, max shift S K, energy moved M erg
     (|.| A) of E erg, line mass m g": the energy the anchor put into the
     lines, i.e. how far the lines drifted from the LES in that step (signed
     and absolute), the total line energy and the line mass of the level;
  - "line vs LES temperature, max |T_line - T_LES| = G K";
  - "write-back: max volume defect V".
The values are averaged over the substeps of each coarse step and plotted
against time. The energy moved is shown per unit mass of the level, as the
mean shift of e over the step [erg/g] and its time integral (the drift the
lines would have accumulated without the anchor).
"""
import argparse
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RE_STEP = re.compile(r"^STEP = (\d+) TIME = ([0-9.eE+-]+) DT")
RE_ANCHOR = re.compile(
    r"CLEM: level (\d) line energy anchored to the LES, max shift ([0-9.eE+-]+) K, "
    r"energy moved ([0-9.eE+-]+) erg \(\|\.\| ([0-9.eE+-]+)\) of ([0-9.eE+-]+) erg, "
    r"line mass ([0-9.eE+-]+) g")
RE_GAP = re.compile(r"CLEM: level (\d) line vs LES temperature, max \|T_line - T_LES\| = ([0-9.eE+-]+)")
RE_VOL = re.compile(r"CLEM: level (\d) write-back: max volume defect ([0-9.eE+-]+)")


def parse(log):
    """Per coarse step: time, and per level the lists of substep values."""
    steps = []
    current = None
    with open(log) as f:
        for line in f:
            m = RE_STEP.match(line)
            if m:
                # the values printed since the previous STEP line belong to
                # the step that ends at this time
                if current is not None:
                    current["t"] = float(m.group(2))
                    steps.append(current)
                current = {"t": None, "lev": {}}
                continue
            if current is None:
                current = {"t": None, "lev": {}}
            for regex, keys in ((RE_ANCHOR, ("shift", "moved", "abs", "line", "mass")),
                                (RE_GAP, ("gap",)), (RE_VOL, ("vol",))):
                m = regex.search(line)
                if m:
                    d = current["lev"].setdefault(int(m.group(1)), {})
                    for k, v in zip(keys, m.groups()[1:]):
                        d.setdefault(k, []).append(float(v))
                    break
    return steps


def series(steps, lev, key, reduce):
    t, y = [], []
    for s in steps:
        d = s["lev"].get(lev, {})
        if key in d:
            t.append(s["t"] * 1.0e6)
            y.append(reduce(d[key]))
    return np.array(t), np.array(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    steps = parse(os.path.join(args.run, "run.log"))
    levels = sorted({lev for s in steps for lev in s["lev"]})
    colours = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c", 3: "#d62728"}

    fig, ax = plt.subplots(5, 1, figsize=(9, 13), sharex=True)
    for lev in levels:
        c = colours.get(lev, "k")
        t, y = series(steps, lev, "shift", max)
        ax[0].plot(t, y, color=c, lw=0.8, label=f"level {lev}")
        # energy moved by the anchor over the step (sum over the substeps),
        # per unit mass of the level: the mean drift of e over the step
        t, moved = series(steps, lev, "moved", sum)
        _, absm = series(steps, lev, "abs", sum)
        _, mass = series(steps, lev, "mass", lambda v: v[-1])
        ax[1].plot(t, moved / mass, color=c, lw=0.8)
        ax[1].plot(t, absm / mass, color=c, lw=0.6, ls=":")
        ax[2].plot(t, np.cumsum(moved / mass), color=c, lw=1.2)
        t, y = series(steps, lev, "gap", max)
        ax[3].plot(t, y, color=c, lw=0.8)
        t, y = series(steps, lev, "vol", max)
        ax[4].plot(t, y, color=c, lw=0.8)

    ax[0].set_ylabel("max anchor shift\nper substep [K]")
    ax[1].set_ylabel("mean e shift per step\n[erg/g] (dotted: |.|)")
    ax[2].set_ylabel("cumulative mean\ne shift [erg/g]")
    ax[3].set_ylabel("max |T_line - T_LES|\nafter anchor [K]")
    ax[4].set_ylabel("max volume defect\nbefore renorm.")
    for a in ax:
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    ax[1].axhline(0.0, color="0.5", lw=0.6)
    ax[2].axhline(0.0, color="0.5", lw=0.6)
    ax[-1].set_xlabel("t [$\\mu$s]")
    ax[0].set_title(f"Line-LES energy drift, {os.path.basename(os.path.normpath(args.run))}, "
                    f"{len(steps)} coarse steps")
    fig.tight_layout()
    out = args.out or os.path.join(args.run, "energy_drift.png")
    fig.savefig(out, dpi=140)
    print(out)


if __name__ == "__main__":
    main()

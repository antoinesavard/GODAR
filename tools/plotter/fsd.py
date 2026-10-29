#######################################
#   FLOE SIZE DISTRIBUTION
# Histograms and cumulative distribution
# of the radii of one or more r*.dat
# files (initial floes or outputs), to
# compare packings.
#######################################
"""
From the repository root:

    python tools/plotter/fsd.py files/r8.dat path/to/rARCTIC.dat \\
        --out plots/plot/packing/fsd.png

Panels: fraction of the floes in each radius bin (bins equally spaced in
log r), fraction of the ice area in each bin, and the cumulative
distribution N(>r) / N with the exponent alpha of a power law
N(>r) ~ r^-alpha fitted between 2 r_min and r_max / 2 (away from the
cut-offs).
"""

import argparse
import pathlib

import numpy as np
import matplotlib.pyplot as plt


def power_law_exponent(r):
    """alpha of N(>r) ~ r^-alpha, least squares in log-log between
    2 r_min and r_max / 2; None if that range is too narrow"""
    r = np.sort(r)[::-1]
    count = np.arange(1, len(r) + 1)
    keep = (r >= 2 * r[-1]) & (r <= r[0] / 2)
    if keep.sum() < 10:
        return None
    slope = np.polyfit(np.log(r[keep]), np.log(count[keep]), 1)[0]
    return -slope


def fsd_figure(radii, labels, nbins=40):
    """figure of the floe size distributions of the lists of radii [m]"""
    radii = [np.asarray(r, dtype=float) / 1e3 for r in radii]   # km
    lo = min(r.min() for r in radii)
    hi = max(r.max() for r in radii)
    bins = np.geomspace(lo, hi * 1.0001, nbins + 1)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    fig.set_layout_engine("tight")
    for r, label in zip(radii, labels):
        number, _ = np.histogram(r, bins)
        area, _ = np.histogram(r, bins, weights=r**2)
        axes[0].stairs(number / len(r), bins, label=label)
        axes[1].stairs(area / area.sum(), bins, label=label)

        s = np.sort(r)[::-1]
        frac = np.arange(1, len(s) + 1) / len(s)
        alpha = power_law_exponent(s)
        line, = axes[2].loglog(
            s, frac,
            label=label + (f", alpha = {alpha:.2f}" if alpha else ""))
        if alpha:
            fit = (s >= 2 * s[-1]) & (s <= s[0] / 2)
            ref = s[fit][len(s[fit]) // 2]
            k = np.searchsorted(-s, -ref)
            axes[2].loglog(s[fit], frac[k] * (s[fit] / ref) ** -alpha, "--",
                           color=line.get_color(), lw=0.8)

    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("radius r [km]")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("fraction of the floes")
    axes[0].set_title("number of floes per bin", fontsize=10)
    axes[1].set_ylabel("fraction of the ice area")
    axes[1].set_title("ice area per bin", fontsize=10)
    axes[2].set_ylabel("N(>r) / N")
    axes[2].set_title("cumulative distribution", fontsize=10)
    for ax in axes:
        ax.legend(fontsize=8)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", type=pathlib.Path,
                        help="files of radii [m], one per line (r*.dat)")
    parser.add_argument("--labels", nargs="+",
                        help="legend of each file (default: name and number "
                        "of floes)")
    parser.add_argument("--out", type=pathlib.Path,
                        help="figure file (default: shown on screen)")
    args = parser.parse_args()

    radii = [np.loadtxt(f, ndmin=1) for f in args.files]
    labels = args.labels or [f"{f.stem} ({len(r)} floes)"
                             for f, r in zip(args.files, radii)]
    fig = fsd_figure(radii, labels)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=200, bbox_inches="tight")
        print(f"wrote {args.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()

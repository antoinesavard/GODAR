#!/usr/bin/env python3
"""
Builds the ARCTIC benchmark case with a given number of floes.

The floes fill the NSIDC sea ice extent of a month (2015-03 by default,
the maximum extent) on the PS-N grid, inside the land mask masks/ps.sdf,
touching each other: their sizes fill the gaps (method fill of
tools/init/initial_packing_position.py), between r_min and 20 r_min as
in CHANNEL, with r_min set so that the --n floes fill the extent (about
91 % of its area is covered, the rest are gaps smaller than r_min, along
the coast and the edge of the extent mostly). In the open areas the
radii are random between 10 and 20 r_min (--seed, the same case on all
machines). The namelist
(tests/namelist/ARCTIC.nml) gets the length and time scales of the
floes: rtree, bond_lim and dt are the CHANNEL values times r_min / 200 m.

The masks and the case are not in the repository: they are built here,
when missing, with the scripts of tools/:

    masks/ps.dat        tools/init/build_masks.py (needs ncdump)
    masks/ps.sdf        tools/init/build_sdf.py (needs numpy and scipy)
    masks/ice_extent/ps_YYYYMM_extent.dat
                        tools/utils/build_ice_extent.py (downloads the
                        extent from NSIDC, needs Pillow)
    the floes           tools/init/initial_packing_position.py (needs
                        numpy, scipy and numba; a few seconds for 5000
                        floes)

From the repository root:

    python tests/bench/make_arctic.py --n 20000 \\
        --out tests/work/bench/ARCTIC_20000
"""

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parents[1]
ROOT = TESTS.parent
MASKS = ROOT / "masks"

sys.path.insert(0, str(TESTS))
import run_tests as rt  # noqa: E402

# CHANNEL values, for floes with r_min = 200 m
R_MIN, RTREE, BOND_LIM, DT = 200.0, 4.5e3, 2e2, 5e-2


def load_tool(path):
    """import a script of tools/ as a module, without running its main"""
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_masks(year, month):
    """PS-N land mask and sea ice extent, built if missing"""
    ps_dat, ps_sdf = MASKS / "ps.dat", MASKS / "ps.sdf"
    if not ps_dat.exists():
        print("building masks/ps.dat")
        tool = load_tool(ROOT / "tools" / "init" / "build_masks.py")
        tool.convert(next(j for j in tool.JOBS if j["out"] == "ps.dat"))
    if not ps_sdf.exists():
        print("building masks/ps.sdf")
        load_tool(ROOT / "tools" / "init" / "build_sdf.py").build_sdf("ps")
    extent = MASKS / "ice_extent" / f"ps_{year:04d}{month:02d}_extent.dat"
    if not extent.exists():
        print(f"building {extent.relative_to(ROOT)}")
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "utils" / "build_ice_extent.py"),
                str(year),
                str(month),
            ],
            check=True,
            cwd=ROOT,
        )
    return ps_sdf, extent


def read_mask(path):
    """0/1 grid of a mask of masks/ (ps.dat, ice extent), row 0 at the
    bottom; returns it with the cell size and the lower left corner [m]"""
    import numpy as np

    with open(path) as f:
        h = f.readline().split()
        nx, ny, dx = int(h[0]), int(h[1]), float(h[2])
        grid = np.fromstring(f.read(), dtype=np.uint8, sep=" ")
    return grid.reshape(ny, nx)[::-1], dx, float(h[3]), float(h[4]) - ny * dx


def check_packing(x, y, r, water, ice, dx, x0, y0):
    """area of the ice extent (water cells only) covered by the floes,
    overlapping pairs of floes and largest overlap (fraction of the
    smaller radius), floes with the centre outside the extent"""
    import numpy as np
    from scipy.spatial import cKDTree

    i, j = (
        cKDTree(np.c_[x, y]).query_pairs(2 * r.max(), output_type="ndarray").T
    )
    gap = np.hypot(x[i] - x[j], y[i] - y[j]) - r[i] - r[j]
    overlap = -gap / np.minimum(r[i], r[j])
    ix = np.clip(((x - x0) // dx).astype(int), 0, ice.shape[1] - 1)
    iy = np.clip(((y - y0) // dx).astype(int), 0, ice.shape[0] - 1)
    return dict(
        coverage=np.pi * np.sum(r**2) / (np.sum(ice & water) * dx**2),
        overlaps=int(np.sum(overlap > 1e-6)),
        max_overlap=overlap.max(initial=0.0),
        outside=int(np.sum(ice[iy, ix] == 0)),
    )


def plot_case(folder, extent):
    """packing.png and packing_fsd.png in the case folder"""
    import matplotlib

    matplotlib.use("Agg")  # no window, under ctest or in a job
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import EllipseCollection
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Rectangle

    x, y, r = (
        np.loadtxt(folder / "files" / f"{v}ARCTIC.dat")
        for v in ("x", "y", "r")
    )
    # ps.dat (1 = water, 0 = land) and the ice extent are on the same grid
    water, dx, x0, y0 = read_mask(MASKS / "ps.dat")
    ice = read_mask(extent)[0]
    stats = check_packing(x, y, r, water, ice, dx, x0, y0)
    summary = (
        f"{len(r)} floes, r = {r.min() / 1e3:.1f} to "
        f"{r.max() / 1e3:.1f} km: {stats['coverage']:.1%} of the ice "
        f"extent covered, {stats['overlaps']} overlapping pairs "
        f"(largest overlap {stats['max_overlap']:.0e} r), "
        f"{stats['outside']} floes with the centre outside the extent"
    )
    print(summary)

    fsd = load_tool(ROOT / "tools" / "plotter" / "fsd.py")
    fig = fsd.fsd_figure([r], [f"{folder.name} ({len(r)} floes)"])
    fig.savefig(folder / "packing_fsd.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # in km
    x, y, r = x / 1e3, y / 1e3, r / 1e3
    ny, nx = water.shape
    xc = (x0 + (np.arange(nx) + 0.5) * dx) / 1e3
    yc = (y0 + (np.arange(ny) + 0.5) * dx) / 1e3
    rows, cols = np.nonzero(ice)
    pad = 2 * r.max()
    view = (
        max(xc[cols.min()] - pad, xc[0]),
        min(xc[cols.max()] + pad, xc[-1]),
        max(yc[rows.min()] - pad, yc[0]),
        min(yc[rows.max()] + pad, yc[-1]),
    )
    k = np.argmin(r)
    half = 2 * r.max()

    fig, (ax_map, ax_zoom) = plt.subplots(
        1,
        2,
        figsize=(14, 7),
        layout="constrained",
        width_ratios=[(view[1] - view[0]) / (view[3] - view[2]), 1],
    )
    for ax, edges in ((ax_map, "none"), (ax_zoom, "k")):
        ax.imshow(
            water,
            origin="lower",
            interpolation="nearest",
            extent=(
                xc[0] - dx / 2e3,
                xc[-1] + dx / 2e3,
                yc[0] - dx / 2e3,
                yc[-1] + dx / 2e3,
            ),
            cmap=ListedColormap(["#d2b48c", "#cfe5f7"]),
            vmin=0,
            vmax=1,
        )
        ax.contour(xc, yc, ice, levels=[0.5], colors="navy", linewidths=0.8)
        ax.add_collection(
            EllipseCollection(
                2 * r,
                2 * r,
                0,
                units="xy",
                offsets=np.c_[x, y],
                offset_transform=ax.transData,
                facecolors="xkcd:crimson",
                edgecolors=edges,
                linewidths=0.3,
                alpha=0.5,
            )
        )
        ax.set_aspect("equal")
        ax.set_xlabel("x [km]")
        ax.set_ylabel("y [km]")
    ax_map.set_xlim(view[0], view[1])
    ax_map.set_ylim(view[2], view[3])
    ax_map.add_patch(
        Rectangle(
            (x[k] - half, y[k] - half),
            2 * half,
            2 * half,
            fill=False,
            edgecolor="k",
        )
    )
    ax_map.set_title(folder.name)
    ax_zoom.set_xlim(x[k] - half, x[k] + half)
    ax_zoom.set_ylim(y[k] - half, y[k] + half)
    ax_zoom.set_title(f"zoom around the smallest floe (r = {r[k]:.1f} km)")
    fig.suptitle(
        f"{summary}\nland: tan, water: blue, ice extent: navy "
        "line, floes: red (darker where they overlap)",
        fontsize=10,
    )
    fig.savefig(folder / "packing.png", dpi=200)
    plt.close(fig)
    print(f"wrote {folder / 'packing.png'} and packing_fsd.png")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--n", type=int, required=True, help="number of floes")
    parser.add_argument(
        "--out", type=Path, required=True, help="folder of the case"
    )
    parser.add_argument(
        "--extent",
        default="2015-03",
        help="month of the NSIDC sea ice extent, YYYY-MM (default: 2015-03)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="seed of the random radii (default: 1)",
    )
    args = parser.parse_args()

    year, month = (int(v) for v in args.extent.split("-"))
    ps_sdf, extent = build_masks(year, month)

    files = args.out.resolve() / "files"
    packed = not all(
        (files / f"{v}ARCTIC.dat").exists()
        for v in ("x", "y", "r", "h", "omega", "theta")
    )
    if packed:
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "tools" / "init" / "initial_packing_position.py"),
                "--mask",
                str(ps_sdf.relative_to(ROOT)),
                "--ice-extent",
                str(extent.relative_to(ROOT)),
                "--n",
                str(args.n),
                "--method",
                "fill",
                "--seed",
                str(args.seed),
                "--expno",
                "ARCTIC",
                "--out",
                str(files),
                "--no-plot",
            ],
            check=True,
            cwd=ROOT,
        )

    mask = files / "maskARCTIC.sdf"
    if not mask.exists():
        mask.symlink_to(ps_sdf)

    radii = [float(v) for v in (files / "rARCTIC.dat").read_text().split()]
    scale = min(radii) / R_MIN
    settings = {
        "rtree": f"{RTREE * scale:.6e}",
        "bond_lim": f"{BOND_LIM * scale:.6e}",
        "dt": f"{DT * scale:.6e}",
    }
    text = rt.set_namelist(
        (TESTS / "namelist" / "ARCTIC.nml").read_text(), settings
    )
    (args.out / "namelist").mkdir(parents=True, exist_ok=True)
    (args.out / "namelist" / "ARCTIC.nml").write_text(text)

    info = dict(
        n=len(radii),
        requested=args.n,
        extent=args.extent,
        seed=args.seed,
        r_min=min(radii),
        r_max=max(radii),
        r_mean=sum(radii) / len(radii),
        scale=scale,
        **{k: float(v) for k, v in settings.items()},
    )
    (args.out / "case.json").write_text(json.dumps(info, indent=4) + "\n")
    print(
        f"ARCTIC case in {args.out}: {len(radii)} floes (of {args.n}), "
        f"r = {min(radii) / 1e3:.1f} to {max(radii) / 1e3:.1f} km, "
        f"rtree = {RTREE * scale / 1e3:.1f} km, dt = {DT * scale:.3g} s"
    )

    if packed or not (files.parent / "packing.png").exists():
        plot_case(files.parent, extent)


if __name__ == "__main__":
    main()

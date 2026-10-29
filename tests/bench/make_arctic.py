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
        --out tests/work/bench/cases/ARCTIC_20000
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
        subprocess.run([sys.executable, str(ROOT / "tools" / "utils"
                                            / "build_ice_extent.py"),
                        str(year), str(month)], check=True, cwd=ROOT)
    return ps_sdf, extent


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                         formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, required=True,
                        help="number of floes")
    parser.add_argument("--out", type=Path, required=True,
                        help="folder of the case")
    parser.add_argument("--extent", default="2015-03",
                        help="month of the NSIDC sea ice extent, YYYY-MM "
                        "(default: 2015-03)")
    parser.add_argument("--seed", type=int, default=1,
                        help="seed of the random radii (default: 1)")
    args = parser.parse_args()

    year, month = (int(v) for v in args.extent.split("-"))
    ps_sdf, extent = build_masks(year, month)

    files = args.out.resolve() / "files"
    if not all((files / f"{v}ARCTIC.dat").exists()
               for v in ("x", "y", "r", "h", "omega", "theta")):
        subprocess.run([sys.executable, str(ROOT / "tools" / "init"
                                            / "initial_packing_position.py"),
                        "--mask", str(ps_sdf.relative_to(ROOT)),
                        "--ice-extent", str(extent.relative_to(ROOT)),
                        "--n", str(args.n), "--method", "fill",
                        "--seed", str(args.seed),
                        "--expno", "ARCTIC", "--out", str(files),
                        "--no-plot"], check=True, cwd=ROOT)

    mask = files / "maskARCTIC.sdf"
    if not mask.exists():
        mask.symlink_to(ps_sdf)

    radii = [float(v) for v in (files / "rARCTIC.dat").read_text().split()]
    scale = min(radii) / R_MIN
    settings = {"rtree": f"{RTREE * scale:.6e}",
                "bond_lim": f"{BOND_LIM * scale:.6e}",
                "dt": f"{DT * scale:.6e}"}
    text = rt.set_namelist((TESTS / "namelist" / "ARCTIC.nml").read_text(),
                           settings)
    (args.out / "namelist").mkdir(parents=True, exist_ok=True)
    (args.out / "namelist" / "ARCTIC.nml").write_text(text)

    info = dict(n=len(radii), requested=args.n, extent=args.extent,
                seed=args.seed, r_min=min(radii), r_max=max(radii),
                r_mean=sum(radii) / len(radii), scale=scale,
                **{k: float(v) for k, v in settings.items()})
    (args.out / "case.json").write_text(json.dumps(info, indent=4) + "\n")
    print(f"ARCTIC case in {args.out}: {len(radii)} floes (of {args.n}), "
          f"r = {min(radii) / 1e3:.1f} to {max(radii) / 1e3:.1f} km, "
          f"rtree = {RTREE * scale / 1e3:.1f} km, dt = {DT * scale:.3g} s")


if __name__ == "__main__":
    main()

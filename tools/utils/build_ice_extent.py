#!/usr/bin/env python3
"""
Build an ice-extent .dat mask from NSIDC Sea Ice Index (G02135).

Run from the repo root:
    python tools/utils/build_ice_extent.py 2015 3

Downloads the monthly extent GeoTIFF from NSIDC, thresholds to ice/no-ice, upsamples 8x to land on the 3.125 km PS-N grid used by ps.dat, and writes:
    masks/ice_extent/ps_{YYYYMM}_extent.dat

Only PS-N is supported here, NSIDC G02135 is native PS-N at 25 km.
"""

import argparse
import pathlib
import sys
import urllib.request

import numpy as np
from PIL import Image

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "masks" / "ice_extent"

MONTH_DIRS = (
    "01_Jan",
    "02_Feb",
    "03_Mar",
    "04_Apr",
    "05_May",
    "06_Jun",
    "07_Jul",
    "08_Aug",
    "09_Sep",
    "10_Oct",
    "11_Nov",
    "12_Dec",
)

# ps.dat target grid
PS = dict(
    nx=2432,
    ny=3584,
    dx=3125.0,
    x_min=-3850000.0,
    y_max=5850000.0,
    projection="PS-N",
)

# NSIDC G02135 native grid (PS-N at 25 km)
NSIDC_NX, NSIDC_NY = 304, 448
UPSAMPLE = PS["nx"] // NSIDC_NX  # = 8


def fetch_geotiff(year: int, month: int) -> pathlib.Path:
    name = f"N_{year:04d}{month:02d}_extent_v4.0.tif"
    local = OUT_DIR / name
    if local.exists():
        print(f"using cached {local}")
        return local
    url = (
        f"https://noaadata.apps.nsidc.org/NOAA/G02135/north/"
        f"monthly/geotiff/{MONTH_DIRS[month - 1]}/{name}"
    )
    print(f"downloading {url}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r:
        local.write_bytes(r.read())
    return local


def read_extent(path: pathlib.Path) -> np.ndarray:
    """NSIDC extent encoding: 0=ocean, 1=ice, 251=pole hole, 253=coast,
    254=missing, 255=land. Treat ice and pole hole as ice."""
    arr = np.array(Image.open(path))
    if arr.shape != (NSIDC_NY, NSIDC_NX):
        sys.exit(
            f"unexpected GeoTIFF shape {arr.shape}, "
            f"expected ({NSIDC_NY}, {NSIDC_NX})"
        )
    return ((arr == 1) | (arr == 251)).astype(np.uint8)


def write_ps_mask(grid: np.ndarray, out_path: pathlib.Path) -> None:
    if grid.shape != (PS["ny"], PS["nx"]):
        sys.exit(f"shape mismatch after upsample: {grid.shape}")
    with out_path.open("w") as f:
        f.write(
            f"{PS['nx']} {PS['ny']} {PS['dx']} "
            f"{PS['x_min']} {PS['y_max']} {PS['projection']}\n"
        )
        np.savetxt(f, grid, fmt="%d")
    print(f"wrote {out_path}  ({PS['nx']} x {PS['ny']})")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("year", type=int, help="e.g. 2015")
    ap.add_argument("month", type=int, choices=range(1, 13))
    args = ap.parse_args()

    tif = fetch_geotiff(args.year, args.month)
    grid = read_extent(tif)
    grid = np.kron(grid, np.ones((UPSAMPLE, UPSAMPLE), dtype=np.uint8))
    out = OUT_DIR / f"ps_{args.year:04d}{args.month:02d}_extent.dat"
    write_ps_mask(grid, out)


if __name__ == "__main__":
    main()

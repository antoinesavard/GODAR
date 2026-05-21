#!/usr/bin/env python3
"""
Build signed-distance-field files from the ASCII mask files.

Run from the repo root:
    python tools/init/build_sdf.py

Reads masks/{ease2,ps}.dat and writes masks/{ease2,ps}.sdf.

SDF convention: positive in water (distance to nearest land in metres),
negative in land (-distance to nearest water in metres). The zero
level-set of the SDF coincides with the mask coast.

Output file format mirrors the .dat masks:
    nx ny dx x_min y_max projection
    <ny rows of nx whitespace-separated floats>

Rows are emitted top-to-bottom, matching the input .dat row order
(y[0] = highest y). The Fortran loader will flip rows so iy=1 sits
at the bottom of the grid.
"""

import pathlib

import numpy as np
from scipy.ndimage import distance_transform_edt

REPO = pathlib.Path(__file__).resolve().parents[2]
MASKS_DIR = REPO / "masks"


def load_mask(path: pathlib.Path):
    with path.open() as f:
        header = f.readline().split()
        nx, ny = int(header[0]), int(header[1])
        dx = float(header[2])
        x_min = float(header[3])
        y_max = float(header[4])
        proj = header[5]
        grid = np.loadtxt(f, dtype=np.uint8).reshape(ny, nx)
    return grid, dx, x_min, y_max, proj


def build_sdf(name: str) -> None:
    src = MASKS_DIR / f"{name}.dat"
    dst = MASKS_DIR / f"{name}.sdf"

    grid, dx, x_min, y_max, proj = load_mask(src)
    ny, nx = grid.shape

    # distance_transform_edt returns Euclidean distance in pixel units
    # from each nonzero cell to the nearest zero cell. Multiplying by dx
    # converts to metres (pixels are isotropic on this mask).
    d_water = (
        distance_transform_edt(grid == 1) * dx
    )  # 0 on land,  > 0 in water
    d_land = distance_transform_edt(grid == 0) * dx  # 0 on water, > 0 on land

    # EDT returns cell-center-to-cell-center distance; the coast actually
    # lies on the cell edge, dx/2 closer. Subtract dx/2 from the magnitude
    # so the zero level-set falls on the water/land cell edge.
    sdf = np.where(
        grid == 1, d_water - dx / 2.0, -(d_land - dx / 2.0)
    ).astype(np.float32)

    with dst.open("w") as f:
        f.write(f"{nx} {ny} {dx} {x_min} {y_max} {proj}\n")
        np.savetxt(f, sdf, fmt="%.2f")

    print(f"wrote {dst}  ({nx} x {ny})")


def main() -> None:
    for name in ("ease2", "ps"):
        build_sdf(name)


if __name__ == "__main__":
    main()

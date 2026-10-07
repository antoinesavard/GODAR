#!/usr/bin/env python3
"""
Crop a rectangular region from masks/ps.dat (NSIDC PS-N 3.125 km) and write
its 0/1 raster + signed-distance field.

Parameters (set below):
    lat_bl, lon_bl : geographic coordinates of the bottom-left corner
                     (degrees). Converted to PS-N (EPSG:3413) via pyproj.
    width, height  : extent of the bounding box in PS-N projection meters.
    region_name    : output basename. Files: masks/{region_name}.dat and
                     masks/{region_name}.sdf.

The bottom-left corner is snapped to the nearest cell edge; width and height
are rounded up to whole cells. The resulting rectangle is therefore at least
as large as requested and aligned with the parent PS-N grid.
"""

import pathlib
import sys

import numpy as np
from pyproj import Transformer
from scipy.ndimage import distance_transform_edt

REPO = pathlib.Path(__file__).resolve().parents[2]
MASKS_DIR = REPO / "masks"

# -------------------------------------
# region parameters
lat_bl = 81.6  # degrees north
lon_bl = -64.7  # degrees east (Greenland Sea side; negative = west)
width = 200.0e3  # PS-N meters along x
height = 200.0e3  # PS-N meters along y
region_name = "nares_strait"
# -------------------------------------

# parent PS-N grid (matches the header of masks/ps.dat written by build_masks.py)
PS_DAT = MASKS_DIR / "ps.dat"
PS_NX = 2432
PS_NY = 3584
PS_DX = 3125.0
PS_X_MIN = -3_850_000.0
PS_Y_MAX = 5_850_000.0
PROJ = "PS-N"


def latlon_to_psn(lat, lon):
    transformer = Transformer.from_crs(
        "EPSG:4326", "EPSG:3413", always_xy=True
    )
    return transformer.transform(lon, lat)


def crop():
    x_bl, y_bl = latlon_to_psn(lat_bl, lon_bl)
    y_top = y_bl + height

    # snap to grid; use floor for the lower edges so the requested point is included
    col_start = int(np.floor((x_bl - PS_X_MIN) / PS_DX))
    row_top = int(np.floor((PS_Y_MAX - y_top) / PS_DX))
    n_cols = int(np.ceil(width / PS_DX))
    n_rows = int(np.ceil(height / PS_DX))

    if col_start < 0 or col_start + n_cols > PS_NX:
        sys.exit(
            f"region extends outside PS-N grid in x: "
            f"cols [{col_start}, {col_start + n_cols}) vs [0, {PS_NX})"
        )
    if row_top < 0 or row_top + n_rows > PS_NY:
        sys.exit(
            f"region extends outside PS-N grid in y: "
            f"rows [{row_top}, {row_top + n_rows}) vs [0, {PS_NY})"
        )

    # snapped origin in PS-N meters
    x_min = PS_X_MIN + col_start * PS_DX
    y_max = PS_Y_MAX - row_top * PS_DX

    print(
        f"requested bottom-left lat/lon = ({lat_bl}, {lon_bl}) "
        f"-> PS-N ({x_bl:.0f}, {y_bl:.0f}) m"
    )
    print(
        f"snapped bounding box: x in [{x_min:.0f}, {x_min + n_cols * PS_DX:.0f}], "
        f"y in [{y_max - n_rows * PS_DX:.0f}, {y_max:.0f}] m"
    )
    print(f"raster: {n_cols} x {n_rows} cells at {PS_DX} m")

    # Direction to the North Pole (which lives at PS-N (0, 0)) from the
    # region centre, expressed two ways:
    #   - bearing: degrees clockwise from +y (i.e., "up" on the plot)
    #              0 = pole directly up, 90 = right, 180 = down, 270 = left
    #   - math angle: degrees counter-clockwise from +x ("right")
    cx = x_min + n_cols * PS_DX / 2.0
    cy = y_max - n_rows * PS_DX / 2.0
    dx_pole = -cx
    dy_pole = -cy
    distance_km = np.hypot(dx_pole, dy_pole) / 1000.0
    bearing = (np.degrees(np.arctan2(dx_pole, dy_pole)) + 360.0) % 360.0
    math_angle = (np.degrees(np.arctan2(dy_pole, dx_pole)) + 360.0) % 360.0
    print(
        f"region centre at PS-N ({cx:.0f}, {cy:.0f}) m, "
        f"{distance_km:.1f} km from the North Pole"
    )
    print(
        f"  to the pole: bearing = {bearing:.1f}° (CW from +y/up), "
        f"math angle = {math_angle:.1f}° (CCW from +x/right)"
    )

    print(f"reading {PS_DAT}...")
    with PS_DAT.open() as f:
        f.readline()  # skip header
        full = np.loadtxt(f, dtype=np.uint8).reshape(PS_NY, PS_NX)

    sub = full[row_top : row_top + n_rows, col_start : col_start + n_cols]
    return sub, n_cols, n_rows, x_min, y_max


def write_ascii(path, grid_top_first, nx, ny, x_min, y_max, fmt):
    with open(path, "w") as f:
        f.write(f"{nx} {ny} {PS_DX} {x_min} {y_max} {PROJ}\n")
        np.savetxt(f, grid_top_first, fmt=fmt)
    print(f"wrote {path}  ({nx} x {ny})")


def main():
    MASKS_DIR.mkdir(parents=True, exist_ok=True)

    water, nx, ny, x_min, y_max = crop()
    write_ascii(
        MASKS_DIR / f"{region_name}.dat", water, nx, ny, x_min, y_max, fmt="%d"
    )

    d_water = distance_transform_edt(water == 1) * PS_DX
    d_land = distance_transform_edt(water == 0) * PS_DX
    # zero level-set on cell edge, not cell center (subtract dx/2)
    sdf = np.where(
        water == 1, d_water - PS_DX / 2.0, -(d_land - PS_DX / 2.0)
    ).astype(np.float32)
    write_ascii(
        MASKS_DIR / f"{region_name}.sdf", sdf, nx, ny, x_min, y_max, fmt="%.2f"
    )


if __name__ == "__main__":
    main()

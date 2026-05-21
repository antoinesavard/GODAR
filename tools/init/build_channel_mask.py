#!/usr/bin/env python3
"""
Build a narrow-channel test mask + SDF.

Geometry (centered on origin):
    - 100 km maximum width, narrowing to a 25 km wide, 100 km long
      channel in the middle. Total y-extent: 200 km.
    - All four grid-edge boundaries are open: water reaches every edge
      of the grid in the wide sections, so particles drifting off any
      edge are simply marked inactive.
    - The channel "shoulders" (the land between |x| = 12.5 km and
      |x| = 50 km at |y| <= 50 km) are the only land cells.

Writes:
    masks/channel.dat   (0/1 raster)
    masks/channel.sdf   (signed-distance field in metres)
"""

import pathlib

import numpy as np
from scipy.ndimage import distance_transform_edt

REPO = pathlib.Path(__file__).resolve().parents[2]
MASKS_DIR = REPO / "masks"

# geometry (metres)
dx = 500.0
half_width_wide = 50.0e3  # max water width / 2 (also half the grid x-extent)
half_width_narrow = 12.5e3  # |x| < this is water inside the channel
channel_y_half = 50.0e3  # channel spans y in [-this, +this]
domain_y_half = 100.0e3  # full grid spans y in [-this, +this]

PROJ = "CHANNEL"


def build():
    Lx = 2.0 * half_width_wide
    Ly = 2.0 * domain_y_half
    nx = int(round(Lx / dx))
    ny = int(round(Ly / dx))

    x_centers = -Lx / 2 + (np.arange(nx) + 0.5) * dx
    y_centers = -Ly / 2 + (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(x_centers, y_centers, indexing="xy")
    aX, aY = np.abs(X), np.abs(Y)

    wide = (aY > channel_y_half) & (aX < half_width_wide)
    narrow = (aY <= channel_y_half) & (aX < half_width_narrow)
    water = (wide | narrow).astype(np.uint8)

    return water, nx, ny, Lx, Ly


def write_ascii(path, grid, nx, ny, Lx, Ly, fmt):
    """Header line + ny rows, top-to-bottom (file row 0 = highest y)."""
    x_min = -Lx / 2
    y_max = Ly / 2
    rows_top_first = grid[::-1]
    with open(path, "w") as f:
        f.write(f"{nx} {ny} {dx} {x_min} {y_max} {PROJ}\n")
        np.savetxt(f, rows_top_first, fmt=fmt)
    print(f"wrote {path}  ({nx} x {ny})")


def main():
    MASKS_DIR.mkdir(parents=True, exist_ok=True)

    water, nx, ny, Lx, Ly = build()
    write_ascii(MASKS_DIR / "channel.dat", water, nx, ny, Lx, Ly, fmt="%d")

    d_water = distance_transform_edt(water == 1) * dx
    d_land = distance_transform_edt(water == 0) * dx
    # zero level-set on the cell edge, not the cell center (subtract dx/2)
    sdf = np.where(
        water == 1, d_water - dx / 2.0, -(d_land - dx / 2.0)
    ).astype(np.float32)
    write_ascii(MASKS_DIR / "channel.sdf", sdf, nx, ny, Lx, Ly, fmt="%.2f")


if __name__ == "__main__":
    main()

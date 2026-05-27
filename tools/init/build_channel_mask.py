#!/usr/bin/env python3
"""
Build a narrow-channel test mask + SDF (asymmetric: channel at the bottom).

Geometry (bottom-anchored; y = 0 is the bottom edge)

Writes:
    masks/channel.dat   (0/1 raster)
    masks/channel.sdf   (signed-distance field in meters)
"""

import pathlib

import numpy as np
from scipy.ndimage import distance_transform_edt

REPO = pathlib.Path(__file__).resolve().parents[2]
MASKS_DIR = REPO / "masks"

# geometry (meters)
dx = 500.0
half_width_wide = 75.0e3  # max water width / 2 (also half the grid x-extent)
half_width_narrow = 12.5e3  # |x| < this is water inside the channel
channel_length = 50.0e3  # channel occupies y in [0, channel_length]
wide_above = 150.0e3  # wide reservoir on top of the channel

PROJ = "CHANNEL"


def build():
    Lx = 2.0 * half_width_wide
    Ly = channel_length + wide_above
    nx = int(round(Lx / dx))
    ny = int(round(Ly / dx))

    x_centers = -Lx / 2 + (np.arange(nx) + 0.5) * dx
    y_centers = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(x_centers, y_centers, indexing="xy")
    aX = np.abs(X)

    channel = (Y <= channel_length) & (aX < half_width_narrow)
    wide = (Y > channel_length) & (aX < half_width_wide)
    water = (channel | wide).astype(np.uint8)

    return water, nx, ny, Lx, Ly


def write_ascii(path, grid, nx, ny, Lx, Ly, fmt):
    """Header line + ny rows, top-to-bottom (file row 0 = highest y)."""
    x_min = -Lx / 2
    y_max = Ly
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
    # zero level-set on cell edge not cell center (subtract dx/2)
    sdf = np.where(
        water == 1, d_water - dx / 2.0, -(d_land - dx / 2.0)
    ).astype(np.float32)
    write_ascii(MASKS_DIR / "channel.sdf", sdf, nx, ny, Lx, Ly, fmt="%.2f")


if __name__ == "__main__":
    main()

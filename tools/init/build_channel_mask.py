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
exp = "12"

# geometry (meters)
dx = 500.0
half_width_wide = 50.0e3  # half the grid x-extent
side_land_width = 0.0e3  # land strip on each side of the wide reservoir (provides confining pressure); set 0 for open sides
half_width_narrow = 12.5e3  # |x| < this is water inside the channel at y=0
channel_length = 50.0e3  # channel occupies y in [0, channel_length]
wide_above = 55.0e3  # wide reservoir on top of the channel
channel_angle_deg = 0.0  # tilt of the shoulder lip from horizontal (degrees). 0 = flat lip at y = channel_length (rectangular shoulders); >0 slopes the lip upward outward, so the shoulders are tallest against the side walls and shortest next to the channel

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

    # Channel column is rectangular (vertical walls); only the shoulder lip
    # tilts. At |x| = half_width_narrow the lip stays at channel_length, and
    # it descends outward by tan(channel_angle_deg) per unit |x|.
    tan_a = np.tan(np.radians(channel_angle_deg))
    res_half = half_width_wide - side_land_width
    shoulder_top = channel_length + np.maximum(0.0, aX - half_width_narrow) * tan_a
    channel = aX < half_width_narrow  # rectangular channel column at any y
    above_lip = (Y >= shoulder_top) & (aX < res_half)
    water = (channel | above_lip).astype(np.uint8)

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
    write_ascii(MASKS_DIR / "{}channel.dat".format(exp), water, nx, ny, Lx, Ly, fmt="%d")

    d_water = distance_transform_edt(water == 1) * dx
    d_land = distance_transform_edt(water == 0) * dx
    # zero level-set on cell edge not cell center (subtract dx/2)
    sdf = np.where(
        water == 1, d_water - dx / 2.0, -(d_land - dx / 2.0)
    ).astype(np.float32)
    write_ascii(MASKS_DIR / "{}channel.sdf".format(exp), sdf, nx, ny, Lx, Ly, fmt="%.2f")


if __name__ == "__main__":
    main()

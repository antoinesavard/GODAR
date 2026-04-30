#!/usr/bin/env python3
"""Plot the ASCII masks produced by tools/init/build_masks.py.

Run from the repo root:
    python tools/plotter/masks.py
Outputs masks/ease2.png and masks/ps.png.
"""
import pathlib

import matplotlib.pyplot as plt
import numpy as np

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
        grid = np.zeros((ny, nx), dtype=np.uint8)
        for j in range(ny):
            grid[j] = np.frombuffer(f.readline().strip().encode(), dtype=np.uint8) - ord("0")
    return grid, dx, x_min, y_max, proj


def plot_mask(name: str):
    grid, dx, x_min, y_max, proj = load_mask(MASKS_DIR / f"{name}.mask")
    ny, nx = grid.shape
    extent = (x_min / 1e3, (x_min + nx * dx) / 1e3,
              (y_max - ny * dx) / 1e3, y_max / 1e3)

    fig, ax = plt.subplots(figsize=(7, 7), dpi=150)
    ax.imshow(grid, origin="upper", extent=extent, cmap="Blues",
              interpolation="nearest", vmin=0, vmax=1)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_title(f"{name}.mask  ({proj}, {nx}x{ny}, dx={dx:.0f} m)")
    ax.set_aspect("equal")

    out = MASKS_DIR / f"{name}.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    for name in ("ease2", "ps"):
        plot_mask(name)

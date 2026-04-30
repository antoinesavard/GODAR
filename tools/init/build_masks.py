#!/usr/bin/env python3
"""Convert NSIDC-0780 sea-ice region NetCDF masks to ASCII 0/1 grids.

Run once from the repo root:
    python tools/init/build_masks.py

Output format (per file):
    nx ny dx x_min y_max projection
    <ny rows of nx chars; '1' = water (sea_ice_region_surface_mask in 0..18),
     '0' = land/ice/fresh-water/disconnected/off-earth/fill>

Rows are emitted top-to-bottom, matching the NetCDF y[0] = highest y.
"""
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
MASKS_DIR = REPO / "masks"

JOBS = [
    dict(nc="NSIDC-0780_SeaIceRegions_EASE2-N3.125km_v1.0.nc",
         out="ease2.mask",
         nx=5760, ny=5760, dx=3125.0,
         x_min=-9_000_000.0, y_max=9_000_000.0,
         projection="EASE2-N"),
    dict(nc="NSIDC-0780_SeaIceRegions_PS-N3.125km_v1.0.nc",
         out="ps.mask",
         nx=2432, ny=3584, dx=3125.0,
         x_min=-3_850_000.0, y_max=5_850_000.0,
         projection="PS-N"),
]


def is_water(v: int) -> bool:
    return 0 <= v <= 18


def convert(job: dict) -> None:
    src = MASKS_DIR / job["nc"]
    dst = MASKS_DIR / job["out"]
    nx, ny = job["nx"], job["ny"]
    total = nx * ny

    proc = subprocess.Popen(
        ["ncdump", "-v", "sea_ice_region_surface_mask", str(src)],
        stdout=subprocess.PIPE, text=True, bufsize=1 << 16,
    )

    try:
        with dst.open("w") as out:
            out.write(
                f'{nx} {ny} {job["dx"]} {job["x_min"]} {job["y_max"]} '
                f'{job["projection"]}\n'
            )
            in_data = False
            row = bytearray()
            cells_done = 0
            for line in proc.stdout:
                if not in_data:
                    eq = line.find("sea_ice_region_surface_mask =")
                    if eq < 0:
                        continue
                    in_data = True
                    line = line[eq:].split("=", 1)[1]
                line = line.strip().rstrip(";").rstrip(",")
                if not line:
                    continue
                for tok in line.split(","):
                    tok = tok.strip()
                    if not tok:
                        continue
                    v = 255 if tok == "_" else int(tok)
                    row.append(0x31 if is_water(v) else 0x30)
                    if len(row) == nx:
                        out.write(row.decode("ascii"))
                        out.write("\n")
                        cells_done += nx
                        del row[:]
                if cells_done >= total:
                    break
    finally:
        proc.stdout.close()
        proc.wait()

    if cells_done != total:
        sys.exit(f"{job['nc']}: parsed {cells_done} cells, expected {total}")
    print(f"wrote {dst}  ({nx} x {ny})")


def main() -> None:
    for j in JOBS:
        convert(j)


if __name__ == "__main__":
    main()

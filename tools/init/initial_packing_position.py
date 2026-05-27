#######################################
#   PACKING
# Recursively places particles either
# in a rectangular domain or inside the
# water region of a signed-distance-field
# mask, following a prescribed floe size
# distribution.
#######################################

import pathlib
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

REPO = pathlib.Path(__file__).resolve().parents[2]
FILES_DIR = REPO / "files"
PLOTS_DIR = REPO / "plots" / "plot" / "packing"

# -------------------------------------
# parameters
domain = (1e5, 5e4)  # fallback (used when mask_file is None)
mask_file = "masks/channel.sdf"
ice_extent_file = None
n = 20000
mu = 800
sigma = 3
r_min = 200
r_max = 4000
dist_name = "pareto"
expno = 11
max_attempts = 3000
small_threshold = 400  # below this, switch to neighbor
# -------------------------------------


# mask / SDF loading
# ------------------------------------------------
def load_sdf(path):
    with open(path) as f:
        h = f.readline().split()
        nx, ny = int(h[0]), int(h[1])
        dx = float(h[2])
        x_min = float(h[3])
        y_max = float(h[4])
        proj = h[5]
        sdf = np.loadtxt(f, dtype=np.float32).reshape(ny, nx)
    # file rows are top-to-bottom (y=y_max first); flip so row 0 = bottom
    sdf = sdf[::-1, :]
    return dict(
        sdf=sdf,
        nx=nx,
        ny=ny,
        dx=dx,
        x_origin=x_min,
        y_origin=y_max - ny * dx,
        Lx=nx * dx,
        Ly=ny * dx,
        proj=proj,
    )


def sdf_at(mask, x, y):
    """Bilinear sample of the SDF at (x, y), treating values as samples
    at cell centers. Returns -inf outside the grid extent so the
    `sdf > r` test naturally rejects off-grid placements."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    scalar = x.ndim == 0
    if scalar:
        x = x[None]
        y = y[None]

    nx, ny, dx = mask["nx"], mask["ny"], mask["dx"]
    inside = (
        (x >= mask["x_origin"])
        & (x <= mask["x_origin"] + nx * dx)
        & (y >= mask["y_origin"])
        & (y <= mask["y_origin"] + ny * dx)
    )
    out = np.full(x.shape, -np.inf, dtype=np.float64)
    if np.any(inside):
        # cell-center frame: fx = 0 at center of column 0
        fx = (x - mask["x_origin"]) / dx - 0.5
        fy = (y - mask["y_origin"]) / dx - 0.5
        ix0 = np.clip(np.floor(fx).astype(int), 0, nx - 2)
        iy0 = np.clip(np.floor(fy).astype(int), 0, ny - 2)
        ix1 = ix0 + 1
        iy1 = iy0 + 1
        tx = np.clip(fx - ix0, 0.0, 1.0)
        ty = np.clip(fy - iy0, 0.0, 1.0)

        s = mask["sdf"]
        ix0c, ix1c = ix0[inside], ix1[inside]
        iy0c, iy1c = iy0[inside], iy1[inside]
        txc = tx[inside]
        tyc = ty[inside]
        out[inside] = (
            (1 - txc) * (1 - tyc) * s[iy0c, ix0c]
            + txc * (1 - tyc) * s[iy0c, ix1c]
            + (1 - txc) * tyc * s[iy1c, ix0c]
            + txc * tyc * s[iy1c, ix1c]
        )
    return out[0] if scalar else out


# spatial hash grid for fast overlap checks
# -------------------------
class SpatialGrid:
    """Bucket the plane into cells of side `cell_size`. A new particle of
    radius r only needs to check the 3x3 cells around its bucket, provided
    cell_size >= 2*r_max."""

    def __init__(self, cell_size):
        self.cell_size = float(cell_size)
        self.cells = {}

    def _key(self, x, y):
        return (int(x // self.cell_size), int(y // self.cell_size))

    def insert(self, x, y, r):
        self.cells.setdefault(self._key(x, y), []).append((x, y, r))

    def has_overlap(self, x, y, r):
        ix, iy = self._key(x, y)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                bucket = self.cells.get((ix + di, iy + dj))
                if not bucket:
                    continue
                for xp, yp, rp in bucket:
                    if (x - xp) ** 2 + (y - yp) ** 2 < (r + rp) ** 2:
                        return True
        return False


# radii distribution
# -------------------------------------------------
def generate_radii(n, dist="normal", params=(0, 1), r_min=1, r_max=10):
    if dist == "lognormal":
        radii = np.random.lognormal(*params, size=n)
    elif dist == "normal":
        radii = np.random.normal(*params, size=n)
        radii = np.where(radii < r_min, r_min, radii)
    elif dist == "uniform":
        radii = np.random.uniform(r_min, r_max, size=n)
    elif dist == "pareto":
        radii = params[0] * np.random.pareto(params[1], size=n)
    else:
        raise ValueError(f"Unsupported distribution: {dist}")
    radii = np.clip(radii, r_min, r_max)
    return np.sort(radii)[::-1]  # largest first


# packing
# ------------------------------------------------------------
def pack_particles(
    radii,
    mask=None,
    ice_mask=None,
    domain_size=None,
    small_threshold=200,
    max_attempts=3000,
):
    if mask is None and domain_size is None:
        raise ValueError("provide either mask or domain_size")

    # bounding box for random sampling: prefer ice_mask, then mask, then domain
    if ice_mask is not None:
        x_lo = ice_mask["x_origin"]
        x_hi = ice_mask["x_origin"] + ice_mask["nx"] * ice_mask["dx"]
        y_lo = ice_mask["y_origin"]
        y_hi = ice_mask["y_origin"] + ice_mask["ny"] * ice_mask["dx"]
    elif mask is not None:
        x_lo, x_hi = mask["x_origin"], mask["x_origin"] + mask["Lx"]
        y_lo, y_hi = mask["y_origin"], mask["y_origin"] + mask["Ly"]
    else:
        x_lo, x_hi = 0.0, domain_size[0]
        y_lo, y_hi = 0.0, domain_size[1]

    grid = SpatialGrid(2.5 * float(np.max(radii)))

    positions = []
    accepted = []
    skipped = 0
    step = max(1, len(radii) // 20)  # progress every ~5%
    t_start = time.time()

    def in_domain(x, y, r):
        if mask is not None and sdf_at(mask, x, y) <= r:
            return False
        if ice_mask is not None and not in_ice_extent(ice_mask, x, y):
            return False
        if mask is None:
            return r <= x <= x_hi - r and r <= y <= y_hi - r
        return True

    for i, r in enumerate(radii):
        if i > 0 and i % step == 0:
            print(
                f"  [{i:>5}/{len(radii)}] placed={len(positions)} "
                f"skipped={skipped} elapsed={time.time() - t_start:.1f}s"
            )
        placed = False
        for _ in range(max_attempts):
            if r > small_threshold or len(positions) < 2:
                x = np.random.uniform(x_lo + r, x_hi - r)
                y = np.random.uniform(y_lo + r, y_hi - r)
            else:
                idx = np.random.choice(len(positions), 2, replace=False)
                p1, p2 = positions[idx[0]], positions[idx[1]]
                r1, r2 = accepted[idx[0]], accepted[idx[1]]
                midpoint = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                d = np.random.normal(size=2)
                d /= np.linalg.norm(d)
                offset = d * ((r1 + r2) / 2 + r) * np.random.uniform(0.9, 1.2)
                x = midpoint[0] + offset[0]
                y = midpoint[1] + offset[1]

            if not in_domain(x, y, r):
                continue
            if grid.has_overlap(x, y, r):
                continue

            positions.append((x, y))
            accepted.append(r)
            grid.insert(x, y, r)
            placed = True
            break

        if not placed:
            skipped += 1

    return np.array(positions), np.array(accepted)


def available_area(mask, ice_mask, domain_size):
    """Area in which a particle can in principle be placed."""
    if ice_mask is not None:
        return int((ice_mask["grid"] == 1).sum()) * ice_mask["dx"] ** 2
    if mask is not None:
        return int((mask["sdf"] > 0).sum()) * mask["dx"] ** 2
    return domain_size[0] * domain_size[1]


def packing_fraction(radii, mask, ice_mask, domain_size):
    covered = float(np.pi * np.sum(radii**2))
    avail = available_area(mask, ice_mask, domain_size)
    return covered / avail if avail > 0 else 0.0


# raster .dat mask loading + helpers
# ------------------------------------------
def load_dat_mask(path):
    """Load an ASCII 0/1 raster mask (same format as ps.dat/ease2.dat).
    The grid is flipped so row 0 sits at the bottom (y increasing)."""
    with open(path) as f:
        h = f.readline().split()
        nx, ny = int(h[0]), int(h[1])
        dx = float(h[2])
        x_min = float(h[3])
        y_max = float(h[4])
        proj = h[5]
        grid = np.loadtxt(f, dtype=np.uint8).reshape(ny, nx)
    return dict(
        grid=grid[::-1, :],
        nx=nx,
        ny=ny,
        dx=dx,
        x_origin=x_min,
        y_origin=y_max - ny * dx,
        proj=proj,
    )


def in_ice_extent(ice_mask, x, y):
    """Cell lookup against ice_mask's own header; outside-grid -> False."""
    ix = int(np.floor((x - ice_mask["x_origin"]) / ice_mask["dx"]))
    iy = int(np.floor((y - ice_mask["y_origin"]) / ice_mask["dx"]))
    if 0 <= ix < ice_mask["nx"] and 0 <= iy < ice_mask["ny"]:
        return ice_mask["grid"][iy, ix] == 1
    return False


def overlay_mask(ax, mask_file_path):
    """Draw the land/water .dat next to the .sdf as the figure background."""
    dat_path = pathlib.Path(str(mask_file_path).replace(".sdf", ".dat"))
    if not dat_path.exists():
        return
    m = load_dat_mask(dat_path)
    extent = (
        m["x_origin"],
        m["x_origin"] + m["nx"] * m["dx"],
        m["y_origin"],
        m["y_origin"] + m["ny"] * m["dx"],
    )
    cmap = ListedColormap(["#d2b48c", "#cfe5f7"])  # land = tan, water = blue
    ax.imshow(
        m["grid"],
        origin="lower",
        extent=extent,
        cmap=cmap,
        interpolation="nearest",
        vmin=0,
        vmax=1,
        zorder=0,
    )
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])


def overlay_ice_extent(ax, ice_mask):
    """Contour the boundary of the ice extent (1 -> 0 transition)."""
    xs = (
        ice_mask["x_origin"]
        + (np.arange(ice_mask["nx"]) + 0.5) * ice_mask["dx"]
    )
    ys = (
        ice_mask["y_origin"]
        + (np.arange(ice_mask["ny"]) + 0.5) * ice_mask["dx"]
    )
    ax.contour(
        xs,
        ys,
        ice_mask["grid"],
        levels=[0.5],
        colors="navy",
        linewidths=1.0,
        zorder=2,
    )


# run
# ----------------------------------------------------------------
t0 = time.time()

if mask_file:
    print(f"loading SDF mask: {mask_file}")
    mask = load_sdf(REPO / mask_file)
    print(
        f"  -> {mask['proj']}, {mask['nx']} x {mask['ny']} cells, "
        f"dx = {mask['dx']:.0f} m, extent "
        f"x=[{mask['x_origin']:.0f}, {mask['x_origin'] + mask['Lx']:.0f}] "
        f"y=[{mask['y_origin']:.0f}, {mask['y_origin'] + mask['Ly']:.0f}]"
    )
else:
    mask = None
    print(
        f"no SDF mask; rectangular domain {domain[0]:.0f} x {domain[1]:.0f} m"
    )

if ice_extent_file:
    print(f"loading ice-extent mask: {ice_extent_file}")
    ice_mask = load_dat_mask(REPO / ice_extent_file)
    cells_in = int(ice_mask["grid"].sum())
    print(
        f"  -> {ice_mask['proj']}, {ice_mask['nx']} x {ice_mask['ny']} cells, "
        f"{cells_in} ice cells ({100*cells_in/ice_mask['grid'].size:.1f}%)"
    )
else:
    ice_mask = None

print(
    f"generating {n} radii ({dist_name}, mu={mu}, sigma={sigma}, "
    f"r in [{r_min}, {r_max}] m)"
)
radii = generate_radii(
    n, dist=dist_name, params=(mu, sigma), r_min=r_min, r_max=r_max
)
print(
    f"  -> min={radii.min():.1f}, max={radii.max():.1f}, "
    f"mean={radii.mean():.1f} m"
)

print(
    f"packing (max_attempts={max_attempts}, small_threshold={small_threshold})..."
)
t_pack = time.time()
positions, accepted_radii = pack_particles(
    radii,
    mask=mask,
    ice_mask=ice_mask,
    domain_size=domain,
    small_threshold=small_threshold,
    max_attempts=max_attempts,
)
print(
    f"packed {len(positions)} / {n} particles in {time.time() - t_pack:.1f} s "
    f"({100*len(positions)/n:.1f}% accepted)"
)

pf = packing_fraction(accepted_radii, mask, ice_mask, domain)
title_extra = f"packing fraction = {pf * 100:.2f}%"
if mask is not None:
    title_extra = f"{mask['proj']}  |  " + title_extra
print(f"packing fraction: {pf * 100:.2f}% of available area")

# figure
# -------------------------------------------------------------
fig, ax = plt.subplots()
fig.set_layout_engine("tight")

if mask is not None:
    overlay_mask(ax, REPO / mask_file)
else:
    ax.set_xlim(0, domain[0])
    ax.set_ylim(0, domain[1])

if ice_mask is not None:
    overlay_ice_extent(ax, ice_mask)

print(f"rendering figure ({len(positions)} circles)...")
for (x, y), r in zip(positions, accepted_radii):
    ax.add_patch(
        plt.Circle((x, y), r, color="xkcd:bright orange", alpha=0.5, zorder=1)
    )

ax.set_aspect("equal")
mean_r_km = float(np.mean(accepted_radii)) / 1e3
params_line = (
    f"expno={expno} | dist={dist_name}(mu={mu}, sigma={sigma}) | "
    f"r in [{r_min/1e3:g}, {r_max/1e3:g}] km | mean r = {mean_r_km:.2f} km"
)
ax.set_title(
    f"{len(positions)} placed / {n} total | {title_extra}\n{params_line}",
    fontsize=9,
)

# size the figure so the axes fill it
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
base_w = 8
fig.set_size_inches(base_w, base_w * abs((ymax - ymin) / (xmax - xmin)))

# save
# ---------------------------------------------------------------
FILES_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

np.savetxt(FILES_DIR / f"r{expno}.dat", accepted_radii, fmt="%16.16f")
np.savetxt(FILES_DIR / f"x{expno}.dat", positions[:, 0], fmt="%16.16f")
np.savetxt(FILES_DIR / f"y{expno}.dat", positions[:, 1], fmt="%16.16f")
np.savetxt(
    FILES_DIR / f"h{expno}.dat", np.ones_like(accepted_radii), fmt="%16.16f"
)
np.savetxt(
    FILES_DIR / f"omega{expno}.dat",
    np.zeros_like(accepted_radii),
    fmt="%16.16f",
)
np.savetxt(
    FILES_DIR / f"theta{expno}.dat",
    np.zeros_like(accepted_radii),
    fmt="%16.16f",
)

fig.savefig(PLOTS_DIR / f"{expno}.png", bbox_inches="tight", dpi=200)
print(f"wrote files/{{x,y,r,h,omega,theta}}{expno}.dat")
print(f"wrote plots/plot/packing/{expno}.png")
print(f"done in {time.time() - t0:.1f} s total")
plt.show()

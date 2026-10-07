import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.transforms import Affine2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cm
import tools.utils.files as tuf
import os
import sys
import pathlib
import sparse

# ----------------------------------------------------------------------
# figures
xaxis_limits = 140  # in km
xoffset = 20
yaxis_limits = 50  # in km
yoffset = 0

trans = False  # transparent background or not
clean = True  # removes the green/red bars
bonds_bool = False  # plots the bonds as rectangles between the disks
bonds_broken = True  # bonds are plotted as red dots in the middle
decay_frames = 1000  # number of frames a red dot stays visible
breaks_min_rel_vel = 0.02  # shows breaks where velocity gradient at break frame exceeds this (km/frame)
breaks_track_midpoint = True  # False: dots frozen at break location
breaks_style = "dot"  # "dot": original midpoint dot, "line": segment between centers, "perp": perpendicular segment at surface midpoint

# possible plots (all mutually exclusive)
bond_num_plot = False  # plots number of bonds per particle
bond_ratio_plot = False  # plots the ratio of fractured bonds per particle, weighted by the size of the particle
thickness = True  # plots thickness fields
stress = False  # plots the stress as facecolor rather than just white
stress_invariant = 1  # J1 or J2 invariant
velocity_x = False  # plots the x-component of velocity
velocity_y = False  # plots the y-component of velocity

# what you want to produce
video = True
image = False

# mask background overlay, overides the axis limits
mask_overlay = True  # draw the mask under the particles
mask_overlay_file = "nares_strait.dat"  # filename inside masks/

# coming from sim
dt = 5e-2  # tstep size in sim
comp = 5e4  # compression in sim

# miscalleneous
output_dir = "../output/NARES_STRAIT/"  # BIG_SIM_20260518
sf = 1e3  # conversion ratio m <-> km
compression = 1  # data compression of videos
start = 0  # starting frame
stop = None  # stopping frame
cbar_horizontal_placement = 0.015  # gap between rightmost subplot and the colorbar (fraction of figure width)

# ----------------------------------------------------------------------

# --------------------------------------
# loading the relevant data
# --------------------------------------

# reading the arguments for the program
try:
    expno = str(sys.argv[1])
    print("expno = {}".format(expno))
except:
    print("No argument provided by sys.")
    expno = str(input("expno = "))
try:
    n = int(sys.argv[2])
    print("number of particles = {}".format(n))
except:
    print("No argument provided by sys.")
    n = int(input("number of particles = "))

print("Reading the files...")

# listing the files to read
filesx = tuf.list_files(output_dir, "x", expno)
filesy = tuf.list_files(output_dir, "y", expno)
filesr = tuf.list_files(output_dir, "r", expno)
filesh = tuf.list_files(output_dir, "h", expno)
filest = tuf.list_files(output_dir, "theta", expno)
fileso = tuf.list_files(output_dir, "omega", expno)
filesb = tuf.list_files(output_dir, "bond", expno)
filestsigxx = tuf.list_files(output_dir, "tsigxx", expno)
filestsigyy = tuf.list_files(output_dir, "tsigyy", expno)
filestsigxy = tuf.list_files(output_dir, "tsigxy", expno)
filestsigyx = tuf.list_files(output_dir, "tsigyx", expno)
filestfx = tuf.list_files(output_dir, "tfx", expno)
filestfy = tuf.list_files(output_dir, "tfy", expno)
filesmom = tuf.list_files(output_dir, "mom", expno)
filesu = tuf.list_files(output_dir, "u", expno) if velocity_x else None
filesv = tuf.list_files(output_dir, "v", expno) if velocity_y else None

# loading the files in memory
x, y, r, h, t, o, b, tfx, tfy, mom = (
    tuf.multiload(output_dir, filesx, 0, n),
    tuf.multiload(output_dir, filesy, 0, n),
    tuf.multiload(output_dir, filesr, 0, n),
    tuf.multiload(output_dir, filesh, 0, n),
    tuf.multiload(output_dir, filest, 0, n),
    tuf.multiload(output_dir, fileso, 0, n),
    tuf.multiload(output_dir, filesb, 1, n),
    tuf.multiload(output_dir, filestfx, 0, n),
    tuf.multiload(output_dir, filestfy, 0, n),
    tuf.multiload(output_dir, filesmom, 0, n),
)
if velocity_x:
    u = tuf.multiload(output_dir, filesu, 0, n)
if velocity_y:
    v = tuf.multiload(output_dir, filesv, 0, n)

# capture full (unsliced) series so break detection runs from frame 0
# regardless of start/stop/compression
if bonds_broken:
    x_full = tuf.check_dim(x / sf)
    y_full = tuf.check_dim(y / sf)
    r_full = tuf.check_dim(r / sf)
    b_full = tuf.check_dim(b, 1)

# compressing the files
x = x[start:stop:compression] / sf
y = y[start:stop:compression] / sf
r = r[start:stop:compression] / sf
h = h[start:stop:compression]
t = t[start:stop:compression]
o = np.sign(o[start:stop:compression])
b = b[start:stop:compression]
tfx = tfx[start:stop:compression]
tfy = tfy[start:stop:compression]
mom = mom[start:stop:compression]
if velocity_x:
    u = u[start:stop:compression] * 100.0  # m/s -> cm/s
if velocity_y:
    v = v[start:stop:compression] * 100.0  # m/s -> cm/s

# check dimensions of the data
x = tuf.check_dim(x)
y = tuf.check_dim(y)
r = tuf.check_dim(r)
h = tuf.check_dim(h)
t = tuf.check_dim(t)
o = tuf.check_dim(o)
b = tuf.check_dim(b, 1)
tfx = tuf.check_dim(tfx)
tfy = tuf.check_dim(tfy)
mom = tuf.check_dim(mom)
if velocity_x:
    u = tuf.check_dim(u)
if velocity_y:
    v = tuf.check_dim(v)

# massaging
t = np.degrees(t)
edge = np.where(o >= 0, "g", "r")
coords = b.coords  # shape (3, nnz)
t_idx = coords[0]
i_idx = coords[1]
j_idx = coords[2]

# --------------------------------------
# precompute broken-bond events (once)
# --------------------------------------
# For each frame k >= 1, find pairs present at k-1 but absent at k, and
# freeze the (px, py) midpoint at the moment of break.
breaks_k = np.empty(0, dtype=np.int32)
breaks_px = np.empty(0)
breaks_py = np.empty(0)
breaks_rel_vel = np.empty(0)
breaks_i = np.empty(0, dtype=np.int64)
breaks_j = np.empty(0, dtype=np.int64)
if bonds_broken:
    print("Precomputing broken-bond events...")
    coords_full = b_full.coords
    order = np.argsort(coords_full[0], kind="stable")
    sc = coords_full[:, order]
    # Break detection runs over [0, stop)
    T_frames = (
        b_full.shape[0] if stop is None else min(int(stop), b_full.shape[0])
    )
    frame_starts = np.searchsorted(sc[0], np.arange(T_frames + 1))

    bk_list, bpx_list, bpy_list, brv_list = [], [], [], []
    bi_list, bj_list = [], []
    prev_set = set(
        zip(
            sc[1, : frame_starts[1]].tolist(),
            sc[2, : frame_starts[1]].tolist(),
        )
    )
    for k in range(1, T_frames):
        s0, s1 = frame_starts[k], frame_starts[k + 1]
        cur_set = set(zip(sc[1, s0:s1].tolist(), sc[2, s0:s1].tolist()))
        broken = prev_set - cur_set
        if broken:
            ij = np.array(list(broken))
            ii, jj = ij[:, 0], ij[:, 1]
            r_i = r_full[k, ii]
            r_j = r_full[k, jj]
            dx = x_full[k, jj] - x_full[k, ii]
            dy = y_full[k, jj] - y_full[k, ii]
            px = x_full[k, ii] + r_i / (r_i + r_j) * dx
            py = y_full[k, ii] + r_i / (r_i + r_j) * dy
            # Velocity gradient
            k_end = T_frames - 1
            duration = k_end - k
            if duration > 0:
                sep_dx_end = x_full[k_end, jj] - x_full[k_end, ii]
                sep_dy_end = y_full[k_end, jj] - y_full[k_end, ii]
                sep_dx_k = dx
                sep_dy_k = dy
                rel_vel = (
                    np.sqrt(
                        (sep_dx_end - sep_dx_k) ** 2
                        + (sep_dy_end - sep_dy_k) ** 2
                    )
                    / duration
                )
            else:
                rel_vel = np.zeros(len(broken))
            bk_list.append(np.full(len(broken), k, dtype=np.int32))
            bpx_list.append(px)
            bpy_list.append(py)
            brv_list.append(rel_vel)
            bi_list.append(ii.astype(np.int64))
            bj_list.append(jj.astype(np.int64))
        prev_set = cur_set

    if bk_list:
        breaks_k = np.concatenate(bk_list)
        breaks_px = np.concatenate(bpx_list)
        breaks_py = np.concatenate(bpy_list)
        breaks_rel_vel = np.concatenate(brv_list)
        breaks_i = np.concatenate(bi_list)
        breaks_j = np.concatenate(bj_list)
    print(
        "  {} bond-break events precomputed across {} frames.".format(
            breaks_k.size, T_frames
        )
    )
    if breaks_rel_vel.size > 0:
        p25, p50, p75, p90, p95, p99 = np.percentile(
            breaks_rel_vel, [25, 50, 75, 90, 95, 99]
        )
        print(
            "  |Δv| percentiles (km/frame): "
            "p25={:.3g} p50={:.3g} p75={:.3g} p90={:.3g} p95={:.3g} p99={:.3g}".format(
                p25, p50, p75, p90, p95, p99
            )
        )
        if breaks_min_rel_vel > 0:
            kept = int(np.sum(breaks_rel_vel >= breaks_min_rel_vel))
            print(
                "  Threshold breaks_min_rel_vel={:g}: keeping {}/{} breaks ({:.1f}%).".format(
                    breaks_min_rel_vel,
                    kept,
                    breaks_k.size,
                    100.0 * kept / breaks_k.size,
                )
            )

# --------------------------------------
# functions for colors
# --------------------------------------


def map_to_color_bond_num(array, cmap=plt.cm.viridis):
    norm = Normalize(vmin=-0.5, vmax=1 * np.nanmax(array) + 0.5)
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    return mapper.to_rgba(array), mapper


def map_to_color_stress(array, cmap=plt.cm.viridis):
    norm = Normalize(vmin=0, vmax=0.01 * np.nanmax(array))
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    return mapper.to_rgba(array), mapper


def map_to_color_bond_ratio(array, cmap=plt.cm.viridis):
    norm = Normalize(vmin=0, vmax=np.nanmax(array))
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    return mapper.to_rgba(array), mapper


def map_to_color_thickness(array, cmap=plt.cm.viridis):
    norm = Normalize(vmin=np.nanmin(array), vmax=np.nanmax(array))
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    return mapper.to_rgba(array), mapper


def map_to_color_velocity_x(array, cmap=plt.cm.coolwarm):
    # symmetric range around zero so 0 maps to the neutral
    norm = Normalize(vmin=-25, vmax=25)
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

    return mapper.to_rgba(array), mapper


def map_to_alpha(array, low, high):
    minimum = np.min(array)
    maximum = np.max(array)

    diff = maximum - minimum
    diffScale = high - low

    return (array - minimum) * (diffScale / diff) + low


alpha = map_to_alpha(h, 0.5, 1)

if bond_num_plot:
    bond_num = (
        sparse.sum(b, axis=1) + sparse.sum(np.transpose(b, (0, 2, 1)), axis=1)
    ).todense()
    cmap = plt.get_cmap(
        "cmo.dense",
        int(np.max(bond_num)) + 1,
    )
    b_cm, mapper = map_to_color_bond_num(bond_num, cmap=cmap)

if bond_ratio_plot:
    bond_num = (
        sparse.sum(b, axis=1) + sparse.sum(np.transpose(b, (0, 2, 1)), axis=1)
    ).todense()
    print(np.sum(bond_num[0]) / 2)
    bond_ratio = np.zeros_like(bond_num)
    # for i in range(bond_num.shape[0]):
    #     bond_ratio[-i - 1] = (
    #         (bond_num[i] - bond_num[-1]) / bond_num[0] / r[-i] * 50 * 100 / sf
    #     )
    for i in range(bond_num.shape[0]):
        bond_ratio[-i - 1] = (
            (bond_num[i] - bond_num[-1]) / r[-i - 1] * np.amin(r[-i - 1])
        )
    cmap = plt.get_cmap("cmo.dense")
    b_cm, mapper = map_to_color_bond_ratio(bond_ratio, cmap=cmap)


if thickness:
    cmap = plt.get_cmap("cmo.dense")
    h = np.where(h > 5, 5, h)
    h_cm, mapper = map_to_color_thickness(h, cmap=cmap)

if velocity_x:
    cmap = plt.get_cmap("coolwarm")
    u_cm, mapper = map_to_color_velocity_x(u, cmap=cmap)

if velocity_y:
    cmap = plt.get_cmap("coolwarm")
    v_cm, mapper = map_to_color_velocity_x(v, cmap=cmap)

if stress:
    # load
    tsigxx, tsigyy, tsigxy, tsigyx = (
        tuf.multiload(output_dir, filestsigxx, 0, n),
        tuf.multiload(output_dir, filestsigyy, 0, n),
        tuf.multiload(output_dir, filestsigxy, 0, n),
        tuf.multiload(output_dir, filestsigyx, 0, n),
    )
    # compress
    tsigxx = tsigxx[start:stop:compression] / sf
    tsigyy = tsigyy[start:stop:compression] / sf
    tsigxy = tsigxy[start:stop:compression] / sf
    tsigyx = tsigyx[start:stop:compression] / sf
    # convert it in real stress
    tsigxx = tsigxx / r**2 / np.pi
    tsigyy = tsigyy / r**2 / np.pi
    tsigxy = tsigxy / r**2 / np.pi
    tsigyx = tsigyx / r**2 / np.pi
    # check dims
    tsigxx = tuf.check_dim(tsigxx)
    tsigyy = tuf.check_dim(tsigyy)
    tsigxy = tuf.check_dim(tsigxy)
    tsigyx = tuf.check_dim(tsigyx)
    # combine
    sigma = np.sqrt(tsigxx**2 + tsigyy**2 + tsigxy**2 + tsigyx**2)
    sigma_cm, mapper = map_to_color_bond_num(sigma)
    dxx = tsigxx - (tsigxx + tsigyy) / 2
    dyy = tsigyy - (tsigxx + tsigyy) / 2
    dxy = tsigxy
    dyx = tsigyx

    j1 = (tsigxx + tsigyy) / 2
    j2 = np.sqrt((dxx**2 + dyy**2 + 2 * dxy**2) / 2)
    j1_cm, mapper = map_to_color_stress(j1)
    j2_cm, mapper = map_to_color_stress(j2)
    j_cm = j1_cm * (2 - stress_invariant) + j2_cm * (stress_invariant - 1)

# --------------------------------------
# compute some things for bonds
# --------------------------------------
if bonds_bool:
    print("Computing the length and orientation of the bonds.")
    lb_data = tuf.lb_func(
        x[t_idx, i_idx], y[t_idx, i_idx], x[t_idx, j_idx], y[t_idx, j_idx]
    )

    angleb_data = tuf.angleb_func(
        x[t_idx, i_idx], y[t_idx, i_idx], x[t_idx, j_idx], y[t_idx, j_idx]
    )

    rb_data = np.minimum(r[t_idx, i_idx], r[t_idx, j_idx])

    lb = sparse.COO(coords, lb_data, shape=b.shape)
    angleb = sparse.COO(coords, angleb_data, shape=b.shape)
    rb = sparse.COO(coords, rb_data, shape=b.shape)
    print("Done")

os.chdir("../plots/anim/")


# --------------------------------------
# some functions for the animation
# --------------------------------------
def init_lists():
    disks = []
    radii = []
    bonds = []
    broken_pairs = []
    num_bonds = np.zeros(n)
    return disks, radii, bonds, broken_pairs, num_bonds


def update_broken_dots(artist, k):
    """
    Update the broken-bond markers for frame k. Style depends on breaks_style:
      "dot"  - one dot at the radius-weighted contact midpoint;
      "line" - one segment between the two particle centers;
      "perp" - one short segment perpendicular to the i->j axis, centered on
               the midpoint between the two surfaces, of length min(r_i, r_j).
    Filtered by the decay window and the breaks_min_rel_vel threshold. When
    breaks_track_midpoint is True the geometry uses live positions/radii at
    frame k; otherwise it is frozen at the break frame.
    """
    is_line = breaks_style in ("line", "perp")

    def clear():
        if is_line:
            artist.set_segments([])
        else:
            artist.set_offsets(np.empty((0, 2)))

    if breaks_k.size == 0:
        clear()
        return
    # k is a sliced display frame; map to the full-series frame index
    # that break detection used. breaks_k holds full-series frame
    # indices.
    full_k = start + k * compression
    keep = (
        (breaks_k > full_k - decay_frames)
        & (breaks_k <= full_k)
        & (breaks_rel_vel >= breaks_min_rel_vel)
    )
    if not np.any(keep):
        clear()
        return
    ii = breaks_i[keep]
    jj = breaks_j[keep]
    kk = breaks_k[keep]

    # positions: live (from full series at full_k) or frozen at the break frame
    if breaks_track_midpoint:
        x_i = x_full[full_k, ii]
        y_i = y_full[full_k, ii]
        x_j = x_full[full_k, jj]
        y_j = y_full[full_k, jj]
    else:
        x_i = x_full[kk, ii]
        y_i = y_full[kk, ii]
        x_j = x_full[kk, jj]
        y_j = y_full[kk, jj]

    # radii: always taken at the break frame so the segment length is fixed
    r_i = r_full[kk, ii]
    r_j = r_full[kk, jj]

    if breaks_style == "dot":
        if breaks_track_midpoint:
            mx = x_i + r_i / (r_i + r_j) * (x_j - x_i)
            my = y_i + r_i / (r_i + r_j) * (y_j - y_i)
        else:
            mx = breaks_px[keep]
            my = breaks_py[keep]
        artist.set_offsets(np.c_[mx, my])
    elif breaks_style in ("line", "perp"):
        dxij = x_j - x_i
        dyij = y_j - y_i
        norm = np.sqrt(dxij * dxij + dyij * dyij)
        norm = np.where(norm == 0.0, 1.0, norm)
        ux_ij = dxij / norm
        uy_ij = dyij / norm

        if breaks_style == "line":
            # frozen length = inter-center distance at the break frame
            seg_len = np.sqrt(
                (x_full[kk, jj] - x_full[kk, ii]) ** 2
                + (y_full[kk, jj] - y_full[kk, ii]) ** 2
            )
            half_len = 0.5 * seg_len
            # midpoint between live centers; segment oriented along live i->j
            mx = 0.5 * (x_i + x_j)
            my = 0.5 * (y_i + y_j)
            dirx = ux_ij
            diry = uy_ij
        else:  # "perp"
            # frozen length = min(r_i, r_j) at the break frame
            half_len = 0.5 * np.minimum(r_i, r_j)
            # midpoint between live surfaces using frozen radii
            mx = 0.5 * (x_i + x_j) + 0.5 * (r_i - r_j) * ux_ij
            my = 0.5 * (y_i + y_j) + 0.5 * (r_i - r_j) * uy_ij
            # segment direction is perpendicular to live i->j
            dirx = -uy_ij
            diry = ux_ij

        e1x = mx - half_len * dirx
        e1y = my - half_len * diry
        e2x = mx + half_len * dirx
        e2y = my + half_len * diry
        segments = np.stack(
            [np.column_stack((e1x, e1y)), np.column_stack((e2x, e2y))],
            axis=1,
        )
        artist.set_segments(segments)

    if decay_frames > 1:
        age = (k - breaks_k[keep]) / max(decay_frames - 1, 1)
        rgba = np.tile(np.array([0.008, 0.0, 0.208, 1.0]), (keep.sum(), 1))
        rgba[:, 3] = np.clip(1.0 - age, 0.0, 1.0)
        if is_line:
            artist.set_colors(rgba)
        else:
            artist.set_facecolors(rgba)


_MASKS_DIR = pathlib.Path(__file__).resolve().parents[2] / "masks"


def apply_mask_overlay(ax):
    if not mask_overlay:
        return
    path = _MASKS_DIR / mask_overlay_file
    with open(path) as f:
        h = f.readline().split()
        nx, ny = int(h[0]), int(h[1])
        dx = float(h[2])
        x_min = float(h[3])
        y_max = float(h[4])
        grid = np.loadtxt(f, dtype=np.uint8).reshape(ny, nx)
    extent = (
        x_min / sf,
        (x_min + nx * dx) / sf,
        (y_max - ny * dx) / sf,
        y_max / sf,
    )
    cmap = ListedColormap(["#d2b48c", "#0a3060"])  # land=tan, water=dark blue
    ax.imshow(
        grid,
        origin="upper",
        extent=extent,
        cmap=cmap,
        interpolation="nearest",
        vmin=0,
        vmax=1,
        zorder=0,
    )
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])


def fit_figure_to_axes(fig, ax, base_width=8, pad=0.02):
    # Grow the figure on each side
    axes = ax if isinstance(ax, (list, tuple)) else [ax]
    xmin, xmax = axes[0].get_xlim()
    ymin, ymax = axes[0].get_ylim()
    aspect = abs((ymax - ymin) / (xmax - xmin))
    # Single axis or vertical stack: scale up height by axis count.
    # Horizontal stack: scale up width by axis count instead.
    if len(axes) <= 1 or stack_vertical:
        fig.set_size_inches(base_width, len(axes) * base_width * aspect)
    else:
        fig.set_size_inches(len(axes) * base_width / aspect, base_width)
    fig.canvas.draw()

    tight_bb = fig.get_tightbbox(fig.canvas.get_renderer())
    fw, fh = fig.get_size_inches()

    # Shrink the figure to *just* the tight bbox + a small pad on each side,
    # rather than preserving matplotlib's default left/bottom gridspec margin.
    new_w = (tight_bb.x1 - tight_bb.x0) + 2 * pad
    new_h = (tight_bb.y1 - tight_bb.y0) + 2 * pad
    shift_x = -tight_bb.x0 + pad
    shift_y = -tight_bb.y0 + pad

    new_positions = []
    for a in fig.axes:
        pos = a.get_position()
        old_x_in = pos.x0 * fw
        old_y_in = pos.y0 * fh
        old_w_in = pos.width * fw
        old_h_in = pos.height * fh
        new_positions.append(
            [
                (old_x_in + shift_x) / new_w,
                (old_y_in + shift_y) / new_h,
                old_w_in / new_w,
                old_h_in / new_h,
            ]
        )

    fig.set_size_inches(new_w, new_h)
    for a, p in zip(fig.axes, new_positions):
        a.set_position(p)
    fig.canvas.draw()

    # Snap auxiliary axes to the position of the right axis
    main_ids = {id(a) for a in axes}
    extras = [a for a in fig.axes if id(a) not in main_ids]
    if extras:
        rightmost = axes[-1].get_position()
        if stack_vertical and len(axes) >= 2:
            topmost = axes[0].get_position()
            cax_y0, cax_h = rightmost.y0, topmost.y1 - rightmost.y0
        else:
            cax_y0, cax_h = rightmost.y0, rightmost.y1 - rightmost.y0
        for a in extras:
            a.set_position(
                [
                    rightmost.x1 + max(0.0, cbar_horizontal_placement),
                    cax_y0,
                    a.get_position().width,
                    cax_h,
                ]
            )
        fig.canvas.draw()


def init_figure(
    trans=False,
    colors=0,
):
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    if trans:
        fig.patch.set_facecolor("None")
    if colors >= 1:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(
            "right",
            size="5%",
            pad=0.1,
        )
    else:
        cax = None

    # ticks
    ax.tick_params(
        which="both",
        direction="out",
        bottom=True,
        top=False,
        left=True,
        right=False,
        labelleft=True,
    )

    return fig, ax, cax


def _displayed_extent():
    """
    Return (width, height) of the data area that will actually be shown.
    When mask_overlay is on, apply_mask_overlay later overwrites xlim/ylim with the mask raster's extent, so use the mask dimensions. Otherwise fall back to the xaxis_limits/yaxis_limits knobs.
    """
    if mask_overlay:
        path = _MASKS_DIR / mask_overlay_file
        with open(path) as f:
            h = f.readline().split()
            nx_m, ny_m = int(h[0]), int(h[1])
        return float(nx_m), float(ny_m)
    return float(xaxis_limits), float(yaxis_limits)


_disp_w, _disp_h = _displayed_extent()
stack_vertical = _disp_w >= _disp_h  # True: 2 rows; False: 2 cols


def init_figure_image(
    trans=False,
    colors=0,
):
    if stack_vertical:
        fig = plt.figure(figsize=(8, 16 * _disp_h / _disp_w), dpi=300)
        gs = fig.add_gridspec(2, 1, width_ratios=[1], wspace=0.05)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    else:
        fig = plt.figure(figsize=(16 * _disp_w / _disp_h, 8), dpi=300)
        gs = fig.add_gridspec(1, 2, height_ratios=[1], wspace=0.05)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)

    # anchor the colorbar to ax1
    if colors >= 1:
        if stack_vertical:
            cax = fig.add_axes(
                [
                    ax1.get_position().x1 + cbar_horizontal_placement,
                    ax1.get_position().y0,
                    0.02,
                    ax0.get_position().y1 - ax1.get_position().y0,
                ]
            )
        else:
            cax = fig.add_axes(
                [
                    ax1.get_position().x1 + cbar_horizontal_placement,
                    ax1.get_position().y0,
                    0.02,
                    ax1.get_position().y1 - ax1.get_position().y0,
                ]
            )
    else:
        cax = None

    ax0.set_aspect("equal")
    ax1.set_aspect("equal")
    if trans:
        fig.patch.set_facecolor("None")

    # ticks
    if stack_vertical:
        ax0.tick_params(
            which="both",
            direction="out",
            bottom=True,
            top=False,
            left=True,
            right=False,
            labelleft=True,
            labelbottom=False,
        )
        ax1.tick_params(
            bottom=True,
            top=False,
            left=True,
            right=False,
            labelleft=True,
        )
    else:
        ax0.tick_params(
            which="both",
            direction="out",
            bottom=True,
            top=False,
            left=True,
            right=False,
            labelleft=True,
            labelbottom=True,
        )
        ax1.tick_params(
            bottom=True,
            top=False,
            left=True,
            right=False,
            labelleft=False,
            labelbottom=True,
        )
    ax = [ax0, ax1]

    return fig, ax, cax


# --------------------------------------------
# video initialization
# --------------------------------------------
if video:
    fig, ax, cax = init_figure(
        trans,
        stress
        + bond_num_plot
        + bond_ratio_plot
        + thickness
        + velocity_x
        + velocity_y,
    )
    ax.set_ylabel(
        r"$y$ [km]",
        rotation=0,
        multialignment="left",
        ha="right",
    )
    ax.set_xlabel(r"$x$ [km]")

    # limits of the plot in kilometers
    ax.set_xlim(-xoffset, xaxis_limits - xoffset)
    ax.set_ylim(0, yaxis_limits)

    # colors
    ax.set_facecolor("xkcd:baby blue")

    apply_mask_overlay(ax)
    if bond_num_plot:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_ticks(
            ticks=np.arange(0, np.max(bond_num) + 1),
            labels=np.arange(0, np.max(bond_num) + 1, dtype=int),
        )
        cb.set_label(
            "Number of\nbonds",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
        )
    elif bond_ratio_plot:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        # cb.set_ticks(
        #     ticks=np.arange(0, np.max(bond_ratio) + 1),
        #     labels=np.arange(0, np.max(bond_ratio) + 1, dtype=int),
        # )
        # cb.set_label(
        #     "fractured\nbonds [%]",
        #     rotation=0,
        #     multialignment="left",
        #     ha="left",
        #     va="top",
        # )
    elif thickness:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "Thickness [m]",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
            position=(0, 0.9),
        )
    elif stress:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "$J_1$ [Pa]",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
        )
    elif velocity_x:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "$u$ [cm/s]",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
            position=(0, 0.9),
        )
    elif velocity_y:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "$v$ [cm/s]",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
            position=(0, 0.9),
        )

    # bonds_broken overlay is independent of the disk-coloring modes above
    if bonds_broken:
        if breaks_style == "dot":
            xrange = ax.get_xlim()[1] - ax.get_xlim()[0]
            scale = 10.0 / xrange
            marker_size = 250.0 * scale**2
            broken_scatter = ax.scatter(
                [],
                [],
                c="xkcd:midnight blue",
                s=marker_size,
                edgecolors="none",
                marker="o",
                zorder=5,
            )
        else:
            broken_scatter = LineCollection(
                [], colors="xkcd:midnight blue", linewidths=1.0, zorder=5
            )
            ax.add_collection(broken_scatter)

    # keep track of time in the figure
    time = fig.text(
        0,
        1.02,
        r"$t = 0\>$hour",
        transform=ax.transAxes,
        horizontalalignment="left",
    )

    # resize the figure so it tightly wraps the axes + decorations
    fit_figure_to_axes(fig, ax)

# --------------------------------------------
# image initialization
# --------------------------------------------
elif image:
    fig, ax, cax = init_figure_image(
        trans,
        stress
        + bond_num_plot
        + bond_ratio_plot
        + thickness
        + velocity_x
        + velocity_y,
    )
    if stack_vertical:
        # both subplots show y label
        ax[0].set_ylabel(
            r"$y$ [km]",
            rotation=0,
            multialignment="left",
            ha="right",
        )
        ax[1].set_ylabel(
            r"$y$ [km]",
            rotation=0,
            multialignment="left",
            ha="right",
        )
        ax[1].set_xlabel(r"$x$ [km]")
    else:
        # both subplots show x label
        ax[0].set_ylabel(
            r"$y$ [km]",
            rotation=0,
            multialignment="left",
            ha="right",
        )
        ax[0].set_xlabel(r"$x$ [km]")
        ax[1].set_xlabel(r"$x$ [km]")

    # limits of the plot in kilometers
    ax[0].set_xlim(-xoffset, xaxis_limits - xoffset)
    ax[0].set_ylim(0, yaxis_limits)
    ax[1].set_xlim(-xoffset, xaxis_limits - xoffset)
    ax[1].set_ylim(0, yaxis_limits)

    # colors
    ax[0].set_facecolor("xkcd:baby blue")
    ax[1].set_facecolor("xkcd:baby blue")

    apply_mask_overlay(ax[0])
    apply_mask_overlay(ax[1])

    if bond_num_plot:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(
            mapper, cax=cax, orientation="vertical", use_gridspec=True
        )
        cb.set_ticks(
            ticks=np.arange(
                0, np.max(bond_num) + 1, (np.max(bond_num) + 1) // 9
            ),
            labels=np.arange(
                0, np.max(bond_num) + 1, (np.max(bond_num) + 1) // 9, dtype=int
            ),
        )
        cb.set_label(
            "Number \nof bonds",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
            position=(0, 0.9),
        )

    elif bond_ratio_plot:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(
            mapper, cax=cax, orientation="vertical", use_gridspec=True
        )
        cb.set_label(
            "Weighted\nfractured\nbonds",
            rotation=0,
            multialignment="left",
            ha="left",
            va="top",
            position=(0, 0.9),
        )

    elif thickness:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(
            mapper, cax=cax, orientation="vertical", extend="max"
        )
        cb.set_label("Thickness [m]", rotation=90)

    elif stress:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "$J_1$ [Pa]",
            rotation=0,
            multialignment="left",
            ha="center",
            va="top",
        )

    elif velocity_x:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label("$u$ [cm/s]", rotation=90)
    elif velocity_y:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label("$v$ [cm/s]", rotation=90)

    # bonds_broken overlay is independent of the disk-coloring modes
    if bonds_broken:
        broken_scatters = []
        for sub_ax in ax:
            if breaks_style == "dot":
                xrange = sub_ax.get_xlim()[1] - sub_ax.get_xlim()[0]
                scale = 10.0 / xrange
                marker_size = 250.0 * scale**2
                sc = sub_ax.scatter(
                    [],
                    [],
                    c="xkcd:midnight blue",
                    s=marker_size,
                    edgecolors="none",
                    marker="o",
                    zorder=5,
                )
            else:
                sc = LineCollection(
                    [], colors="xkcd:midnight blue", linewidths=1.0, zorder=5
                )
                sub_ax.add_collection(sc)
            broken_scatters.append(sc)

    # keep track of time in the figure
    time0 = fig.text(
        0,
        1.02,
        r"$t = 0\>$hour",
        transform=ax[0].transAxes,
        horizontalalignment="left",
    )
    time1 = fig.text(
        0,
        1.02,
        r"$t = 0\>$hour",
        transform=ax[1].transAxes,
        horizontalalignment="left",
    )

    # resize the figure so it tightly wraps the axes + decorations
    fit_figure_to_axes(fig, ax)

# --------------------------------------
# functions for the animation/imagination
# --------------------------------------
disks, radii, bonds, broken_pairs, num_bonds = init_lists()


def init(ax, time):
    print("Initial drawing in process")
    for i in range(x.shape[-1]):
        p = np.array([x[0, i], y[0, i]])
        disk, rad = tuf.draw(ax, p, r[0, i], t[0, i], edge[0, i])
        if i == len(x[-1]) - 1:
            disks.append(disk)
            radii.append(rad)
            continue
        disks.append(disk)
        radii.append(rad)
        if (
            bond_num_plot
            or bond_ratio_plot
            or thickness
            or stress
            or velocity_x
            or velocity_y
        ):
            continue
        if bonds_bool:
            for j in range(i + 1, b.shape[-1]):
                if b[0, j, i]:
                    bond = tuf.draw_bond(
                        ax,
                        p,
                        lb[0, j, i],
                        angleb[0, j, i],
                        radius=2 * rb[0, j, i],
                    )
                    bonds.append(bond)
                    num_bonds[i] += 1
    time.set_text("")
    return disks, radii, bonds


def animate(k, time):
    if k % 50 == 0:
        print("Frame: {}".format(k))
    mask_k = b.coords[0] == k
    i_bonds = b.coords[1, mask_k]
    j_bonds = b.coords[2, mask_k]
    for i, (disk, rad) in enumerate(zip(disks, radii)):
        p = np.array([x[k, i], y[k, i]])
        disk.center = p
        disk.radius = r[k, i]
        disk.set_alpha(1)
        if bond_num_plot or bond_ratio_plot:
            disk.set_facecolor(b_cm[k, i])
            disk.set_linewidth(0)
        if thickness:
            disk.set_facecolor(h_cm[k, i])
            disk.set_linewidth(0)
        if stress:
            disk.set_facecolor(j_cm[k, i])
            disk.set_linewidth(0)
        if velocity_x:
            disk.set_facecolor(u_cm[k, i])
            disk.set_linewidth(0)
        if velocity_y:
            disk.set_facecolor(v_cm[k, i])
            disk.set_linewidth(0)
        rad.xy = p
        rad.angle = t[k, i]
        rad.set_width(r[k, i])
        rad.set_edgecolor(edge[k, i])
        if clean is True:
            rad.set_visible(False)
        # if i == len(disks) - 1:
        #     continue
        if (
            bond_num_plot
            or bond_ratio_plot
            or thickness
            or stress
            or velocity_x
            or velocity_y
        ):
            continue
        if bonds_bool:
            for n, bond in enumerate(bonds):
                i = i_bonds[n]
                j = j_bonds[n]

                if b[k, i, j] is False:
                    bond.set_visible(False)
                else:
                    theta = np.deg2rad(angleb[k, i, j])
                    width = lb[k, i, j]
                    height = 2 * rb[k, i, j]
                    # unit vector from i -> j
                    ux = np.cos(theta)
                    uy = np.sin(theta)

                    # contact point
                    px = x[k, i] + height * uy / 2
                    py = y[k, i] - height * ux / 2

                    bond.angle = angleb[k, i, j]
                    bond.set_width(width)
                    bond.set_height(height)
                    bond.xy = (px, py)

    # broken-bond dots
    if bonds_broken and not (bond_num_plot or bond_ratio_plot):
        update_broken_dots(broken_scatter, k)

    time.set_text(
        r"$t = {:.0f}\>$hour".format(
            dt * comp * (start + k * compression) / 3600.0
        )
    )

    if k == 0 or k == r.shape[0] - 1:
        fig.savefig(
            "../../plots/plot/collision{}-{}.png".format(expno, k + 1),
            bbox_inches="tight",
            pad_inches=0.02,
        )
        fig.savefig(
            "../../plots/plot/collision{}-{}.pdf".format(expno, k + 1),
            bbox_inches="tight",
            pad_inches=0.02,
        )
        fig.savefig(
            "../../plots/plot/collision{}-{}_trans.png".format(expno, k + 1),
            transparent=True,
            bbox_inches="tight",
            pad_inches=0.02,
        )
        fig.savefig(
            "../../plots/plot/collision{}-{}_trans.pdf".format(expno, k + 1),
            transparent=True,
            bbox_inches="tight",
            pad_inches=0.02,
        )

    return disks, radii, bonds


# --------------------------------------
# function for plots
# --------------------------------------
def imaginate(
    k,
    time,
    ax,
    scatter=None,
):
    print("Initial drawing in process")
    for i in range(x.shape[-1]):
        p = np.array([x[k, i], y[k, i]])
        disk, rad = tuf.draw(ax, p, r[k, i], t[k, i], edge[k, i])

        if bond_num_plot or bond_ratio_plot:
            disk.set_facecolor(b_cm[k, i])
            disk.set_linewidth(0)

        if thickness:
            disk.set_facecolor(h_cm[k, i])
            disk.set_linewidth(0)

        if stress:
            disk.set_facecolor(j_cm[k, i])
            disk.set_linewidth(0)

        if velocity_x:
            disk.set_facecolor(u_cm[k, i])
            disk.set_linewidth(0)

        if velocity_y:
            disk.set_facecolor(v_cm[k, i])
            disk.set_linewidth(0)

        if clean is True:
            rad.set_visible(False)

        if i == len(x[-1]) - 1:
            disks.append(disk)
            radii.append(rad)
            continue

        if (
            bond_num_plot
            or bond_ratio_plot
            or thickness
            or stress
            or velocity_x
            or velocity_y
        ):
            continue

        if bonds_bool:
            for j in range(0, b.shape[-1]):
                if b[k, j, i]:
                    bond = tuf.draw_bond(
                        ax,
                        p,
                        lb[k, j, i],
                        angleb[k, j, i],
                        radius=2 * rb[k, j, i],
                    )
                    bonds.append(bond)
                    num_bonds[i] += 1
                    if bond_num_plot or bond_ratio_plot:
                        bond.set_visible(False)
        disks.append(disk)
        radii.append(rad)

    # broken-bond dots
    if bonds_broken and not (bond_num_plot or bond_ratio_plot):
        update_broken_dots(scatter, k)

    time.set_text(
        r"$t = {:.0f}\>$hour".format(
            dt * comp * (start + k * compression) / 3600.0
        )
    )

    return disks, radii, bonds


# --------------------------------------
# the wrappers for the animation
# --------------------------------------
# define some wrappers
def init_wrapper():
    return init(ax, time)


def animate_wrapper(k):
    return animate(k, time)


# --------------------------------------
# animating and saving the data
# --------------------------------------

if video:
    anim = FuncAnimation(
        fig,
        animate_wrapper,
        frames=x.shape[0],
        init_func=init_wrapper,
        interval=10,
        repeat=False,
        blit=False,
    )
    print("Animating the disks for your eyes.")
    tuf.save_or_show_animation(
        anim, 1, "../../plots/anim/collision{}.mp4".format(expno)
    )

if image:
    imaginate(
        0, time0, ax[0], scatter=broken_scatters[0] if bonds_broken else None
    )
    imaginate(
        x.shape[0] - 1,
        time1,
        ax[1],
        scatter=broken_scatters[1] if bonds_broken else None,
    )
    fig.savefig(
        "../../plots/plot/collision{}-startstop.png".format(expno),
        bbox_inches="tight",
        pad_inches=0.02,
    )
    fig.savefig(
        "../../plots/plot/collision{}-startstop.pdf".format(expno),
        bbox_inches="tight",
        pad_inches=0.02,
    )
    fig.savefig(
        "../../plots/plot/collision{}-startstop_trans.png".format(expno),
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.02,
    )
    fig.savefig(
        "../../plots/plot/collision{}-startstop_trans.pdf".format(expno),
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.02,
    )

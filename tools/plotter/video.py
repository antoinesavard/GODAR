import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cm
import tools.utils.files as tuf
import os
import sys
import sparse

# ----------------------------------------------------------------------
# figures
xaxis_limits = 5  # in km
xoffset = 0
yaxis_limits = 2  # in km
yoffset = 0
trans = True  # transparent background or not
clean = False  # removes the green/red bars
bonds_bool = True  # plots the bonds as rectangles between the disks
bonds_broken = False  # bonds are plotted as red dots in the middle

# possible plots (all mutually exclusive)
bond_num_plot = False  # plots number of bonds per particle
bond_ratio_plot = False  # plots the ratio of fractured bonds per particle, weighted by the size of the particle
thickness = False  # plots thickness fields
stress = False  # plots the stress as facecolor rather than just white
stress_invariant = 2  # J1 or J2 invariant

# what you want to produce
video = True
image = False

# coming from sim
dt = 1e-3  # tstep size in sim
comp = 1e5  # compression in sim

# miscalleneous
output_dir = "../output/"
sf = 1e3  # conversion ratio m <-> km
compression = 1  # data compression of videos
start = 0  # starting frame
stop = None  # stopping frame
cbar_horizontal_placement = (
    -0.01
)  # horizontal placement of the colorbar in the image, in fraction of the axis width

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

# massaging
t = np.degrees(t)
edge = np.where(o >= 0, "g", "r")
coords = b.coords  # shape (3, nnz)
t_idx = coords[0]
i_idx = coords[1]
j_idx = coords[2]

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
    h = np.where(h > 2, 2, h)
    h_cm, mapper = map_to_color_thickness(h, cmap=cmap)

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

    j1 = (dxx + dyy) / 2
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


def init_figure(
    trans=False,
    colors=0,
):
    fig, ax = plt.subplots()
    fig.set_layout_engine("tight")
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


def init_figure_image(
    trans=False,
    colors=0,
):
    fig = plt.figure(figsize=(8, 16 * yaxis_limits / xaxis_limits), dpi=300)
    gs = fig.add_gridspec(2, 1, width_ratios=[1], wspace=0.05)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)

    # anchor the colorbar to ax1
    if colors >= 1:
        cax = fig.add_axes(
            [
                ax1.get_position().x1 + cbar_horizontal_placement,
                ax1.get_position().y0,
                0.02,
                ax0.get_position().y1 - ax1.get_position().y0,
            ]
        )
    else:
        cax = None

    ax0.set_aspect("equal")
    ax1.set_aspect("equal")
    if trans:
        fig.patch.set_facecolor("None")

    # ticks
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
    ax = [ax0, ax1]

    return fig, ax, cax


# --------------------------------------------
# video initialization
# --------------------------------------------
if video:
    fig, ax, cax = init_figure(
        trans,
        stress + bond_num_plot + bond_ratio_plot + thickness,
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
    if bond_num_plot:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_ticks(
            ticks=np.arange(0, np.max(bond_num) + 1),
            labels=np.arange(0, np.max(bond_num) + 1, dtype=int),
        )
        cb.set_label(
            "Number of\nbonds", rotation=0, multialignment="left", ha="left", va="top"
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
            "$J_1$ [Pa]", rotation=0, multialignment="left", ha="left", va="top"
        )
    elif bonds_broken:
        xrange = ax.get_xlim()[1] - ax.get_xlim()[0]
        ref_range = 10
        ref_size = 15
        scale = ref_range / xrange
        marker_size = ref_size * scale**2
        broken_scatter = ax.scatter(
            [], [], c="xkcd:bright orange", s=marker_size, zorder=5
        )

    # keep track of time in the figure
    time = fig.text(0, 1.02, "", transform=ax.transAxes, horizontalalignment="left")

# --------------------------------------------
# image initialization
# --------------------------------------------
elif image:
    fig, ax, cax = init_figure_image(
        trans,
        stress + bond_num_plot + bond_ratio_plot + thickness,
    )
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

    # limits of the plot in kilometers
    ax[0].set_xlim(-xoffset, xaxis_limits - xoffset)
    ax[0].set_ylim(0, yaxis_limits)
    ax[1].set_xlim(-xoffset, xaxis_limits - xoffset)
    ax[1].set_ylim(0, yaxis_limits)

    # colors
    ax[0].set_facecolor("xkcd:baby blue")
    ax[1].set_facecolor("xkcd:baby blue")

    if bond_num_plot:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical", use_gridspec=True)
        cb.set_ticks(
            ticks=np.arange(0, np.max(bond_num) + 1, (np.max(bond_num) + 1) // 9),
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
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical", use_gridspec=True)
        # cb.set_ticks(
        #     ticks=np.arange(0, np.max(bond_ratio) + 1, (np.max(bond_ratio) + 1)),
        #     labels=np.arange(
        #         0, np.max(bond_ratio) + 1, (np.max(bond_ratio) + 1), dtype=int
        #     ),
        # )
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
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical", extend="max")
        cb.set_label(
            "Thickness [m]",
            rotation=90,
            multialignment="left",
            ha="center",
            va="top",
            position=(0, 0.5),
        )

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

    elif bonds_broken:
        xrange = ax[1].get_xlim()[1] - ax[1].get_xlim()[0]
        ref_range = 10
        ref_size = 15
        scale = ref_range / xrange
        marker_size = ref_size * scale**2
        broken_scatter = ax[1].scatter(
            [], [], c="xkcd:bright orange", s=marker_size, zorder=5
        )

    # keep track of time in the figure
    time0 = fig.text(0, 1.02, "", transform=ax[0].transAxes, horizontalalignment="left")
    time1 = fig.text(0, 1.02, "", transform=ax[1].transAxes, horizontalalignment="left")

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
        if bond_num_plot or bond_ratio_plot or thickness or stress:
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
        rad.xy = p
        rad.angle = t[k, i]
        rad.set_width(r[k, i])
        rad.set_edgecolor(edge[k, i])
        if clean is True:
            rad.set_visible(False)
        # if i == len(disks) - 1:
        #     continue
        if bond_num_plot or bond_ratio_plot or thickness or stress:
            continue
        if bonds_bool:
            for n, bond in enumerate(bonds):
                i = i_bonds[n]
                j = j_bonds[n]

                if b[k, i, j] is False:
                    bond.set_visible(False)
                else:
                    pb = np.array(
                        [
                            x[k, i] + r[k, i] * np.sin(np.deg2rad(angleb[k, i, j])),
                            y[k, i] - r[k, i] * np.cos(np.deg2rad(angleb[k, i, j])),
                        ]
                    )
                    bond.xy = pb
                    bond.angle = angleb[k, i, j] * b[k, i, j]
                    bond.set_width(lb[k, i, j] * b[k, i, j])
                    bond.set_height(2 * rb[k, i, j] * b[k, i, j])

        elif bonds_broken:
            # Detect bond breaking
            if k > 0:
                mask_k = b.coords[0] == k
                i_k = b.coords[1, mask_k]
                j_k = b.coords[2, mask_k]

                mask_prev = b.coords[0] == 0
                i_prev = b.coords[1, mask_prev]
                j_prev = b.coords[2, mask_prev]

                prev_bonds = set(zip(i_prev, j_prev))
                current_bonds = set(zip(i_k, j_k))
                broken = prev_bonds - current_bonds

                if broken:
                    i_idx = np.array([p[0] for p in broken])
                    j_idx = np.array([p[1] for p in broken])

                    # radii of the two particles
                    r_i = r[k, i_idx]
                    r_j = r[k, j_idx]

                    # vector from i -> j
                    dx = x[k, j_idx] - x[k, i_idx]
                    dy = y[k, j_idx] - y[k, i_idx]

                    # weighted midpoint along the bond
                    px = x[k, i_idx] + r_i / (r_i + r_j) * dx
                    py = y[k, i_idx] + r_i / (r_i + r_j) * dy

                    broken_scatter.set_offsets(np.c_[px, py])
                else:
                    broken_scatter.set_offsets(np.empty((0, 2)))

    time.set_text(
        r"$t = {}\>$hour".format(round(dt * comp * compression * (k + 1) / 60 / 60))
    )

    if k == 0 or k == r.shape[0] - 1:
        fig.savefig("../../plots/plot/collision{}-{}.png".format(expno, k + 1))
        fig.savefig("../../plots/plot/collision{}-{}.pdf".format(expno, k + 1))

    return disks, radii, bonds


# --------------------------------------
# function for plots
# --------------------------------------
def imaginate(
    k,
    time,
    ax,
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

        if clean is True:
            rad.set_visible(False)

        if i == len(x[-1]) - 1:
            disks.append(disk)
            radii.append(rad)
            continue

        if bond_num_plot or bond_ratio_plot or thickness or stress:
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
        elif bonds_broken:
            # Detect bond breaking
            if k > 0:
                mask_k = b.coords[0] == k
                i_k = b.coords[1, mask_k]
                j_k = b.coords[2, mask_k]

                mask_prev = b.coords[0] == 0
                i_prev = b.coords[1, mask_prev]
                j_prev = b.coords[2, mask_prev]

                prev_bonds = set(zip(i_prev, j_prev))
                current_bonds = set(zip(i_k, j_k))
                broken = prev_bonds - current_bonds

                if broken:
                    i_idx = np.array([p[0] for p in broken])
                    j_idx = np.array([p[1] for p in broken])

                    # radii of the two particles
                    r_i = r[k, i_idx]
                    r_j = r[k, j_idx]

                    # vector from i -> j
                    dx = x[k, j_idx] - x[k, i_idx]
                    dy = y[k, j_idx] - y[k, i_idx]

                    # weighted midpoint along the bond
                    px = x[k, i_idx] + r_i / (r_i + r_j) * dx
                    py = y[k, i_idx] + r_i / (r_i + r_j) * dy

                    broken_scatter.set_offsets(np.c_[px, py])
                else:
                    broken_scatter.set_offsets(np.empty((0, 2)))

        disks.append(disk)
        radii.append(rad)

    time.set_text(
        r"$t = {}\>$hour".format(round(dt * comp * compression * (k + 1) / 60 / 60))
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
    imaginate(0, time0, ax[0])
    imaginate(x.shape[0] - 1, time1, ax[1])
    fig.savefig("../../plots/plot/collision{}-startstop.png".format(expno))
    fig.savefig("../../plots/plot/collision{}-startstop.pdf".format(expno))

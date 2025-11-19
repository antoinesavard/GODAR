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

# ----------------------------------------------------------------------
# figures
xaxis_limits = 140  # in km
xoffset = 20
yaxis_limits = 50  # in km
yoffset = 0
trans = True  # transparent background or not
clean = True  # removes the green/red bars

# possible plots
bond_num_plot = False  # plots number of bonds per particle
thickness = False  # plots thickness fields
stress = False  # plots the stress as facecolor rather than just white
stress_invariant = 2  # J1 or J2 invariant

# what you want
video = True
image = False

# coming from sim
dt = 1e-3  # tstep size in sim
comp = 1e5  # compression in sim

# miscalleneous
sf = 1e3  # conversion ratio m <-> km
compression = 1  # data compression of videos

# ----------------------------------------------------------------------

# --------------------------------------
# loading the relevant data
# --------------------------------------

# reading the arguments for the program
output_dir = "../output/"
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

# loading the files in memory
x, y, r, h, t, o, b = (
    tuf.multiload(output_dir, filesx, 0, n),
    tuf.multiload(output_dir, filesy, 0, n),
    tuf.multiload(output_dir, filesr, 0, n),
    tuf.multiload(output_dir, filesh, 0, n),
    tuf.multiload(output_dir, filest, 0, n),
    tuf.multiload(output_dir, fileso, 0, n),
    tuf.multiload(output_dir, filesb, 1, n),
)

# compressing the files
x = x[::compression] / sf
y = y[::compression] / sf
r = r[::compression] / sf
h = h[::compression]
t = t[::compression]
o = np.sign(o[::compression])
b = b[::compression]

# check dimensions of the data
x = tuf.check_dim(x)
y = tuf.check_dim(y)
r = tuf.check_dim(r)
h = tuf.check_dim(h)
t = tuf.check_dim(t)
o = tuf.check_dim(o)
b = tuf.check_dim(b, 1)

# massaging
t = np.degrees(t)
edge = np.where(o >= 0, "g", "r")
lb = np.zeros_like(b)
rb = np.zeros_like(b)
angleb = np.zeros_like(b)

# --------------------------------------
# functions for colors
# --------------------------------------


def map_to_color(array, cmap=plt.cm.viridis):
    norm = Normalize(vmin=-0.5, vmax=1 * np.nanmax(array) + 0.5)
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
    bond_num = np.sum(b, axis=1) + np.sum(np.transpose(b, (0, 2, 1)), axis=1)
    bond_ratio = np.zeros_like(bond_num)
    for i in range(bond_num.shape[0]):
        bond_ratio[-i - 1] = bond_num[i] - bond_num[-1]
    cmap = plt.get_cmap(
        "cmo.dense",
        int(np.max(bond_num)) + 1,
    )
    b_cm, mapper = map_to_color(bond_num, cmap=cmap)

if thickness:
    cmap = plt.get_cmap("cmo.dense")
    h_cm, mapper = map_to_color(h, cmap=cmap)

if stress:
    # load
    tsigxx, tsigyy, tsigxy, tsigyx = (
        tuf.multiload(output_dir, filestsigxx, 0, n),
        tuf.multiload(output_dir, filestsigyy, 0, n),
        tuf.multiload(output_dir, filestsigxy, 0, n),
        tuf.multiload(output_dir, filestsigyx, 0, n),
    )
    # compress
    tsigxx = tsigxx[::compression] / sf
    tsigyy = tsigyy[::compression] / sf
    tsigxy = tsigxy[::compression] / sf
    tsigyx = tsigyx[::compression] / sf
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
    sigma_cm, mapper = map_to_color(sigma)
    dxx = tsigxx - (tsigxx + tsigyy) / 2
    dyy = tsigyy - (tsigxx + tsigyy) / 2
    dxy = tsigxy
    dyx = tsigyx

    j1 = (dxx + dyy) / 2
    j2 = np.sqrt((dxx**2 + dyy**2 + 2 * dxy**2) / 2)
    j1_cm, mapper = map_to_color(j1)
    j2_cm, mapper = map_to_color(j2)
    j_cm = j1_cm * (2 - stress_invariant) + j2_cm * (stress_invariant - 1)

# --------------------------------------
# compute some things for bonds
# --------------------------------------
if not clean:
    print("Compute the length and orientation of the bonds...")
    for i in range(b.shape[-1] - 1):
        for j in range(i + 1, b.shape[-2]):
            lb[:, i, j] = tuf.lb_func(x[:, i], y[:, i], x[:, j], y[:, j])
            angleb[:, i, j] = tuf.angleb_func(x[:, i], y[:, i], x[:, j], y[:, j])
    for i in range(rb.shape[0] - 1):
        rb[i] = tuf.rb_func(r[i], r[i])

os.chdir("../plots/anim/")

# --------------------------------------
# create the figure to animate
# --------------------------------------


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
    # fig, ax = plt.subplots(1, 2, sharey=True)
    fig = plt.figure(layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)

    # anchor the colorbar to ax1
    if colors >= 1:
        cax = inset_axes(
            ax1,
            width="100%",  # colorbar width
            height="5%",  # match ax1 height
            loc="lower center",
            bbox_to_anchor=(-0.55, -0.3, 1, 1),
            bbox_transform=ax1.transAxes,
            borderpad=0,
        )
    else:
        cax = None
    # fig = plt.figure("constrained")
    # if colors >= 1:
    #     gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.03], wspace=0.05)
    #     ax0 = fig.add_subplot(gs[0, 0])
    #     ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)
    #     cax = fig.add_subplot(gs[0, 2])
    # else:
    #     gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
    #     ax0 = fig.add_subplot(gs[0, 0])
    #     ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)
    #     cax = None

    # fig.set_layout_engine("tight")
    ax0.set_aspect("equal")
    ax1.set_aspect("equal")
    if trans:
        fig.patch.set_facecolor("None")
    # if colors >= 1:
    #     divider = make_axes_locatable(ax[1])
    #     cax = divider.append_axes(
    #         "right",
    #         size="5%",
    #         pad=0.1,
    #     )
    # else:
    #     cax = None

    # ticks
    ax0.tick_params(
        which="both",
        direction="out",
        bottom=True,
        top=False,
        left=True,
        right=False,
        labelleft=True,
    )
    ax1.tick_params(
        bottom=True,
        top=False,
        left=False,
        right=False,
        labelleft=False,
    )
    ax = [ax0, ax1]

    return fig, ax, cax


if video:
    fig, ax, cax = init_figure(
        trans,
        stress + bond_num_plot + thickness,
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
    elif thickness:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "Thickness [m]", rotation=0, multialignment="left", ha="left", va="top"
        )
    elif stress:
        ax.set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="vertical")
        cb.set_label(
            "$J_1$ [Pa]", rotation=0, multialignment="left", ha="left", va="top"
        )
    # keep track of time in the figure
    time = fig.text(0, 1.02, "", transform=ax.transAxes, horizontalalignment="left")
elif image:
    fig, ax, cax = init_figure_image(
        trans,
        stress + bond_num_plot + thickness,
    )
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
    if bond_num_plot:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="horizontal")
        cb.set_ticks(
            ticks=np.arange(0, np.max(bond_num) + 1),
            labels=np.arange(0, np.max(bond_num) + 1, dtype=int),
        )
        cb.set_label(
            "Number of bonds", rotation=0, multialignment="left", ha="center", va="top"
        )
    elif thickness:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="horizontal")
        cb.set_label(
            "Thickness [m]", rotation=0, multialignment="left", ha="center", va="top"
        )
    if stress:
        ax[0].set_facecolor("white")
        ax[1].set_facecolor("white")
        cb = fig.colorbar(mapper, cax=cax, orientation="horizontal")
        cb.set_label(
            "$J_1$ [Pa]",
            rotation=0,
            multialignment="left",
            ha="center",
            va="top",
        )
    # keep track of time in the figure
    time0 = fig.text(0, 1.02, "", transform=ax[0].transAxes, horizontalalignment="left")
    time1 = fig.text(0, 1.02, "", transform=ax[1].transAxes, horizontalalignment="left")

# second figure
fig_strip = plt.figure(dpi=300, figsize=(4 * xaxis_limits / yaxis_limits, 4))
ax_strip = fig_strip.add_axes([0, 0, 1, 1])

# limits of the plot in kilometers
ax_strip.set_xlim(0, xaxis_limits)
ax_strip.set_ylim(0, yaxis_limits)

# keep track of time in the figure
time_strip = fig_strip.text(2, 2, "", transform=ax_strip.transAxes)


def init_lists():
    disks = []
    radii = []
    bonds = []
    num_bonds = np.zeros(n)
    return disks, radii, bonds, num_bonds


# --------------------------------------
# the functions for the animation
# --------------------------------------

disks, radii, bonds, num_bonds = init_lists()


def init(ax, time):
    print("Initial drawing in process")
    for i in range(x.shape[-1]):
        p = np.array([x[0, i], y[0, i]])
        disk, rad = tuf.draw(ax, p, r[0, i], t[0, i], edge[0, i])
        if i == len(x[-1]) - 1:
            disks.append(disk)
            radii.append(rad)
            continue
        for j in range(i + 1, b.shape[-1]):
            if b[0, i, j]:
                bond = tuf.draw_bond(
                    ax,
                    p,
                    lb[0, i, j],
                    angleb[0, i, j],
                    radius=2 * rb[0, i, j],
                )
                bonds.append(bond)
                num_bonds[i] += 1
                if bond_num_plot:
                    bond.set_visible(False)
        disks.append(disk)
        radii.append(rad)
    time.set_text("")
    return disks, radii, bonds


def animate(k, time):
    if k % 50 == 0:
        print("Frame: {}".format(k))
    loc = np.argwhere(b[0] > 0)
    # loc = np.argwhere((b[k] + b[k + 1] == 2) | (b[k] - b[k + 1] == -1))
    # lost = np.argwhere(b[k] - b[k + 1] == 1)
    # gain = np.argwhere(b[k] - b[k + 1] == -1)
    # if (lost.shape[0] >= 1) | (gain.shape[0] >= 1):
    #     num_bonds[lost[:, 0]] -= 1
    #     num_bonds[gain[:, 0]] += 1
    for i, (disk, rad) in enumerate(zip(disks, radii)):
        p = np.array([x[k, i], y[k, i]])
        disk.center = p
        disk.radius = r[k, i]
        disk.set_alpha(1)
        if bond_num_plot:
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
        if i == len(disks) - 1:
            continue
        if bond_num_plot:
            continue
        index = np.sum(num_bonds[:i], dtype=int)
        for j, bond in enumerate(
            bonds
        ):  # need to change this otherwise I loop over all bonds, and I should only be looping on the bonds attached to disk i, there should be a counter that keeps track of where we are in the list
            if i == 0:
                index = 0
            if j < index or j >= index + num_bonds[i]:
                continue
            elif b[k, loc[j, 0], loc[j, 1]] is False:
                bond.set_visible(False)
            else:
                pb = np.array(
                    [
                        x[k, i]
                        + r[k, i] * np.sin(np.deg2rad(angleb[k, loc[j, 0], loc[j, 1]])),
                        y[k, i]
                        - r[k, i] * np.cos(np.deg2rad(angleb[k, loc[j, 0], loc[j, 1]])),
                    ]
                )
                bond.xy = pb
                bond.angle = (
                    angleb[k, loc[j, 0], loc[j, 1]] * b[k, loc[j, 0], loc[j, 1]]
                )
                bond.set_width(lb[k, loc[j, 0], loc[j, 1]] * b[k, loc[j, 0], loc[j, 1]])
                bond.set_height(
                    2 * rb[k, loc[j, 0], loc[j, 1]] * b[k, loc[j, 0], loc[j, 1]]
                )
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
        if bond_num_plot:
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
        if bond_num_plot:
            continue
        for j in range(i + 1, b.shape[-1]):
            if b[k, i, j]:
                bond = tuf.draw_bond(
                    ax,
                    p,
                    lb[k, i, j],
                    angleb[k, i, j],
                    radius=2 * rb[k, i, j],
                )
                bonds.append(bond)
                num_bonds[i] += 1
                if bond_num_plot:
                    bond.set_visible(False)
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


def init_wrapper_strip():
    return init(ax_strip, time_strip)


def animate_wrapper(k):
    return animate(k, time)


def animate_wrapper_strip(k):
    return animate(k, time_strip)


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

# anim_strip = FuncAnimation(
#     fig_strip,
#     animate_wrapper_strip,
#     frames=x.shape[0],
#     init_func=init_wrapper_strip,
#     interval=10,
#     repeat=False,
#     blit=False,
# )

# print("Animating the disks for the computer's eyes.")
# tuf.save_or_show_animation(
#     anim_strip, 1, "../../plots/anim/strip-collision{}.mp4".format(expno)
# )

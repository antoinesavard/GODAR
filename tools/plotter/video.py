import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tools.utils.files as tuf
import os
import sys

# ----------------------------------------------------------------------
# figures
xaxis_limits = 6  # in km
xoffset = 0
yaxis_limits = 3  # in km
xlenght = 6.2  # 2.4
trans = False  # transparent background or not

# coming from sim
nt = 2e5
dt = 1e-3  # tstep size in sim
comp = 1e2  # compression in sim

# miscalleneous
sf = 1e3  # conversion ratio m <-> km
compression = 5  # data compression of videos

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
h = h[::compression] / sf
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
angleb = np.zeros_like(b)

# --------------------------------------
# compute some things
# --------------------------------------
print("Compute the length and orientation of the bonds...")
for i in range(b.shape[-1] - 1):
    for j in range(i + 1, b.shape[-2]):
        lb[:, i, j] = tuf.lb_func(x[:, i], y[:, i], x[:, j], y[:, j])
        angleb[:, i, j] = tuf.angleb_func(x[:, i], y[:, i], x[:, j], y[:, j])

os.chdir("../plots/anim/")

# --------------------------------------
# create the figure to animate
# --------------------------------------


def init_figure(number_plots, xaxis, yaxis, x=2.15, y=3.2, trans=False):
    # init plot
    basefigsizex = 7
    figsizex = x * (0.325581 + 0.474419 + 0.1 * (number_plots + 1))
    figsizey = y
    fig = plt.figure(
        dpi=300,
        figsize=(figsizex, figsizey),
    )
    if trans:
        fig.patch.set_facecolor("None")

    # definitions for the axis

    # graphs
    extraspaceinit = 0.1 * basefigsizex / figsizex
    graphsizex = (0.9 - extraspaceinit) / number_plots
    graphsizey = yaxis / xaxis * graphsizex * figsizex / figsizey
    separation = (1 - extraspaceinit - number_plots * graphsizex) / (number_plots + 1)

    bottom, height = (0.125, graphsizey)
    left, width = (separation + extraspaceinit, graphsizex)

    rect = [
        left,
        bottom,
        width,
        height,
    ]

    ax = fig.add_axes(rect)

    # ticks
    ax.tick_params(
        which="both",
        direction="out",
        bottom=False,
        top=False,
        left=True,
        right=False,
        labelleft=True,
    )

    return fig, ax


fig, ax = init_figure(
    1,
    xaxis_limits,
    yaxis_limits,
    xlenght,
    3.2,
    trans,
)

fig.text(0.03, 0.5, r"$y$ [km]")
fig.text(0.54, 0.025, r"$x$ [km]")

# limits of the plot in kilometers
ax.set_xlim(-xoffset, xaxis_limits - xoffset)
ax.set_ylim(0, yaxis_limits)

# colors
ax.set_facecolor("xkcd:baby blue")

# second figure
fig_strip = plt.figure(dpi=300, figsize=(4 * xaxis_limits / yaxis_limits, 4))
ax_strip = fig_strip.add_axes([0, 0, 1, 1])

# limits of the plot in kilometers
ax_strip.set_xlim(0, xaxis_limits)
ax_strip.set_ylim(0, yaxis_limits)

# keep track of time in the figure
time = fig.text(0.65, 1.02, "", transform=ax.transAxes)
time_strip = fig_strip.text(2, 2, "", transform=ax_strip.transAxes)

disks = []
radii = []
bonds = []
num_bonds = np.zeros(n)

# --------------------------------------
# the functions for the animation
# --------------------------------------


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
                    radius=2 * min(r[0, i], r[0, j]),
                )
                bonds.append(bond)
                num_bonds[i] += 1
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
        rad.xy = p
        rad.angle = t[k, i]
        rad.set_width(r[k, i])
        rad.set_edgecolor(edge[k, i])
        if i == len(disks) - 1:
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
                bond.set_height(2 * min(r[k, i], r[k, j]))
    time.set_text("t = {}s".format(round(dt * comp * compression * (k + 1))))

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
tuf.save_or_show_animation(anim, 1, "../../plots/anim/collision{}.mp4".format(expno))

anim_strip = FuncAnimation(
    fig_strip,
    animate_wrapper_strip,
    frames=x.shape[0],
    init_func=init_wrapper_strip,
    interval=10,
    repeat=False,
    blit=False,
)

print("Animating the disks for the computer's eyes.")
tuf.save_or_show_animation(
    anim_strip, 1, "../../plots/anim/strip-collision{}.mp4".format(expno)
)

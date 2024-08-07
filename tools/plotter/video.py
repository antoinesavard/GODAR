import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import tools.utils.files as ff
import os
import sys

# ----------------------------------------------------------------------

xaxis_limits = 30  # in km
yaxis_limits = 30  # in km
sf = 1e3  # conversion ratio m <-> km
compression = 1  # data compression

# ----------------------------------------------------------------------

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
filesx = ff.list_files(output_dir, "x", expno)
filesy = ff.list_files(output_dir, "y", expno)
filesr = ff.list_files(output_dir, "r", expno)
filesh = ff.list_files(output_dir, "h", expno)
filest = ff.list_files(output_dir, "theta", expno)
fileso = ff.list_files(output_dir, "omega", expno)
filesb = ff.list_files(output_dir, "bond", expno)

# loading the files in memory
x, y, r, h, t, o, b = (
    ff.multiload(output_dir, filesx),
    ff.multiload(output_dir, filesy),
    ff.multiload(output_dir, filesr),
    ff.multiload(output_dir, filesh),
    ff.multiload(output_dir, filest),
    ff.multiload(output_dir, fileso),
    ff.multiload(output_dir, filesb, 1, n),
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
x = ff.check_dim(x)
y = ff.check_dim(y)
r = ff.check_dim(r)
h = ff.check_dim(h)
t = ff.check_dim(t)
o = ff.check_dim(o)
b = ff.check_dim(b, 1)

# massaging
t = np.degrees(t)
edge = np.where(o >= 0, "g", "r")
lb = np.zeros_like(b)
angleb = np.zeros_like(b)

print("Compute the length and orientation of the bonds...")
for i in range(b.shape[-1] - 1):
    for j in range(i + 1, b.shape[-1]):
        lb[:, i, j] = ff.lb_func(x[:, i], y[:, i], x[:, j], y[:, j])
        angleb[:, i, j] = ff.angleb_func(x[:, i], y[:, i], x[:, j], y[:, j])

os.chdir("../plots/anim/")

fig = plt.figure(dpi=300, figsize=(4 * xaxis_limits / yaxis_limits, 4))
ax = fig.add_axes([0.14, 0.14, 0.8, 0.8])
ax.set_ylabel("Position [km]", rotation=90)
ax.set_xlabel("Position [km]")

# limits of the plot in kilometers
ax.set_xlim(0, xaxis_limits)
ax.set_ylim(0, yaxis_limits)

# second figure
fig_strip = plt.figure(dpi=300, figsize=(4 * xaxis_limits / yaxis_limits, 4))
ax_strip = fig_strip.add_axes([0, 0, 1, 1])

# limits of the plot in kilometers
ax_strip.set_xlim(0, xaxis_limits)
ax_strip.set_ylim(0, yaxis_limits)

# keep track of time in the figure
time = ax.text(0.8, 1.02, "", transform=ax.transAxes)
time_strip = ax_strip.text(2, 2, "", transform=ax_strip.transAxes)

disks = []
radii = []
bonds = []
num_bonds = np.zeros(n)


def init(ax, time):
    print("Initial drawing in process")
    for i in range(x.shape[-1]):
        p = np.array([x[0, i], y[0, i]])
        disk, rad = ff.draw(ax, p, r[0, i], t[0, i], edge[0, i])
        if i == len(x[-1]) - 1:
            disks.append(disk)
            radii.append(rad)
            continue
        for j in range(i + 1, b.shape[-1]):
            if b[0, i, j]:
                bond = ff.draw_bond(
                    ax,
                    p,
                    lb[0, i, j],
                    angleb[0, i, j],
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
                bond.xy = p
                bond.angle = (
                    angleb[k, loc[j, 0], loc[j, 1]] * b[k, loc[j, 0], loc[j, 1]]
                )
                bond.set_width(lb[k, loc[j, 0], loc[j, 1]] * b[k, loc[j, 0], loc[j, 1]])
    time.set_text("t = {}".format(k + 1))

    return disks, radii, bonds


# define some wrappers
def init_wrapper():
    return init(ax, time)


def init_wrapper_strip():
    return init(ax_strip, time_strip)


def animate_wrapper(k):
    return animate(k, time)


def animate_wrapper_strip(k):
    return animate(k, time_strip)


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
ff.save_or_show_animation(anim, 1, "../../plots/anim/collision{}.mp4".format(expno))

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
ff.save_or_show_animation(
    anim_strip, 1, "../../plots/anim/strip-collision{}.mp4".format(expno)
)

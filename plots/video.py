import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from files import *
import os
import sys

output_dir = "output/"
try:
    expno = str(sys.argv[1])
except:
    print("No argument provided by sys.")
    expno = str(input("expno = "))
sf = 1e3
compression = 1

filesx = list_files(output_dir, "x", expno)
filesy = list_files(output_dir, "y", expno)
filesr = list_files(output_dir, "r", expno)
filesh = list_files(output_dir, "h", expno)
filest = list_files(output_dir, "theta", expno)
fileso = list_files(output_dir, "omega", expno)

x, y, r, h, t, o = (
    multiload(output_dir, filesx),
    multiload(output_dir, filesy),
    multiload(output_dir, filesr),
    multiload(output_dir, filesh),
    multiload(output_dir, filest),
    multiload(output_dir, fileso),
)

x = x[::compression]
y = y[::compression]
r = r[::compression]
h = h[::compression]
t = t[::compression]
o = np.sign(o[::compression])
edge = np.where(o >= 0, "g", "r")

os.chdir("plots/")

fig = plt.figure(dpi=300, figsize=(4, 4))
ax = fig.add_axes([0.14, 0.14, 0.8, 0.8])

time = ax.text(0.02, 0.95, "", transform=ax.transAxes)

ax.set_xlim(0, 10 * sf)
ax.set_ylim(0, 10 * sf)

disks = []
radii = []


def init():
    for i in range(len(r[-1])):
        p = np.array([x[0, i], y[0, i]])
        disk, rad = draw(ax, p, r[0, i], t[0, i], edge[0, i])
        disks.append(disk)
        radii.append(rad)
    time.set_text("")
    return disks


def animate(i):
    if i % 50 == 0:
        print("Frame: {}".format(i))
    for j, (disk, rad) in enumerate(zip(disks, radii)):
        p = np.array([x[i, j], y[i, j]])
        disk.center = p
        rad.xy = p
        rad.angle = 0  # t[i, j]
        rad.edgecolor = edge[i, j]
    time.set_text("iteration = {}".format(i))
    return disks


anim = FuncAnimation(
    fig,
    animate,
    frames=x.shape[0],
    init_func=init,
    interval=10,
    repeat=False,
    blit=False,
)

save_or_show_animation(anim, 1, "collision{}.mp4".format(expno))

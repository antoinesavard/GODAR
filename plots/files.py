import numpy as np
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def list_files(directory: str, datatype: str, expno: str) -> list:
    """
    Gives all files with given extension inside given directory.

    Args:
        directory (str): directory name
        extension (str): file extension to list

    Returns:
        (list): all files with given extension inside given directory
    """
    from os import path, listdir

    files = (f for f in listdir(directory) if f.startswith(datatype) and f.endswith(expno))

    return [f for f in files]


def multiload(output_dir, files: list) -> np.ndarray:

    data = []
    for file in files:
        fic = open(output_dir + file)
        data.append(np.loadtxt(fic))
        fic.close()

    return np.stack(data, axis=0)


def draw(ax, r, radius):
    """Add this Particle's Circle patch to the Matplotlib Axes ax."""

    circle = Circle(xy=r, radius=radius, edgecolor="b", fill=True)
    ax.add_patch(circle)
    return circle


def save_or_show_animation(anim, save, filename="collision.mp4"):
    if save:
        Writer = animation.writers["ffmpeg"]
        writer = Writer(fps=60, bitrate=1800)
        anim.save(filename, writer=writer)
    else:
        plt.show()

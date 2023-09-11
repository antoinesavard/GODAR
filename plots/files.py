import numpy as np
from matplotlib.patches import Circle, Rectangle
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

    files = (
        f
        for f in listdir(directory)
        if f.startswith(datatype) and f.endswith(expno)
    )

    return [f for f in files]


def multiload(output_dir, files: list, bond=0, n=None) -> np.ndarray:
    if not bond:
        data = []
        for file in files:
            fic = open(output_dir + file)
            data.append(np.loadtxt(fic))
            fic.close()

        data = np.stack(data, axis=0)

        return data[0] if data.shape[0] == 1 else data

    elif bond:
        print("Reading bonds...")
        for file in files:
            with open(output_dir + file) as fic:
                data = np.loadtxt(fic).reshape(-1, n, n)

        return data[0] if data.shape[0] == 1 else data


def draw(ax, r, radius, angle, edge):
    """Add this Particle's Circle patch to the Matplotlib Axes ax."""
    circle = Circle(xy=r, radius=radius, edgecolor="b", fill=True, zorder=0)
    angle_rect = Rectangle(
        xy=r,
        height=1 / 1000,
        width=radius,
        angle=angle,
        edgecolor=edge,
        fill=True,
    )
    ax.add_patch(circle)
    ax.add_patch(angle_rect)
    return circle, angle_rect


def draw_bond(ax, r, lb, angleb):
    bond_rect = Rectangle(
        xy=r,
        height=1 / 1000,
        width=lb,
        angle=angleb,
        edgecolor=(0, 0, 0, 0.5),
        fill=True,
        zorder=2,
    )
    ax.add_patch(bond_rect)
    return bond_rect


def save_or_show_animation(anim, save, filename="collision.mp4"):
    if save:
        Writer = animation.writers["ffmpeg"]
        writer = Writer(fps=60, bitrate=-1)
        anim.save(filename, writer=writer)
    else:
        plt.show()


def lb_func(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def angleb_func(x1, y1, x2, y2):
    return np.degrees(np.arctan2(y2 - y1, x2 - x1))

import numpy as np
import matplotlib.pyplot as plt

def void_ratio_plot(vr_means_list: list, vr_stds_list: list, expno: str):
    """
    Creates a time series of the void ratio evolution

    Args:
    - vr_means_list (list): list of the mean of void ratios evolving in time.
    - vr_stds_list (list): list of standard deviation of void ratios evolving in time.
    - expno (str): experience number for when we are saving the figures.
    """

    fig = plt.figure(dpi=300, figsize=(6, 4))

    # definitions for the axes
    left, width = 0.14, 0.8
    bottom, height = 0.12, 0.8
    rect_scatter = [left, bottom, width, height]

    ax = fig.add_axes(rect_scatter)

    # plot of the mean
    ax.plot(vr_means_list[1:])

    # Adding shaded region for standard deviation
    ax.fill_between(np.arange(vr_means_list[1:].shape[0]), vr_means_list[1:] - 3 * vr_stds_list[1:], vr_means_list[1:] + 3 * vr_stds_list[1:], color='blue', alpha=0.2)
    
    # axis labels
    ax.set_xlabel("Time step number")
    ax.set_ylabel("Void ratio")

    fig.savefig("../plots/plot/{}void_ratio.pdf".format(expno))
    fig.savefig("../plots/plot/{}void_ratio.png".format(expno))

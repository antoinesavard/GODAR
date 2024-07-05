# ------------------------------------------------
# import modules
# ------------------------------------------------

import numpy as np
import cv2
from numba import njit, prange

# ------------------------------------------------
# define some global variables
# ------------------------------------------------

num_iteration = 100

# ------------------------------------------------
# define the functions
# ------------------------------------------------


# Function to check if a point is in a disk
@njit
def in_disk(x, y, image):
    x, y = int(x), int(y)
    color = image[y, x]
    lower_blue = np.array([0, 0, 0])
    upper_blue = np.array([200, 200, 255])
    white = np.array([245, 245, 245])
    return (
        np.all(lower_blue <= color)
        and np.all(color <= upper_blue)
        and not np.all(color > white)
    )


def find_first_uniform_column(matrix):
    for col in range(matrix.shape[1]):
        column_data = matrix[:, col]
        first_value = column_data[0]
        total = sum([first_value == item for item in column_data])
        if sum(total) == 3:
            return col
    return -1  # Return -1 if no such column exists


@njit(parallel=True)
def monte_carlo(num_samples, frame_rgb, frame_width, frame_height, frame_max):
    """
    Function that generates points randomly in the domain and check whether they are inside disks or not.

    Args:
        num_samples (_type_): _description_
        frame_rgb (_type_): _description_
        frame_width (_type_): _description_
        frame_height (_type_): _description_
        frame_max (_type_): _description_

    Returns:
        _type_: _description_
    """

    points_in_disks = 0
    for _ in prange(num_samples):
        x = np.random.uniform(0, frame_width - (frame_width - frame_max))
        y = np.random.uniform(0, frame_height)
        if in_disk(x, y, frame_rgb):
            points_in_disks += 1

    return points_in_disks


def outer_loop(num_samples, frame_rgb, frame_width, frame_height, frame_max):
    # function that loops over each frame in order to get significant statistics
    void_ratio_arr = np.zeros(num_iteration)

    for i in range(num_iteration):

        # # here
        # fig = plt.figure(dpi=100, figsize=(8, 8))
        # # definitions for the axes
        # left, width = 0, 1
        # bottom, height = 0, 1
        # rect_scatter = [left, bottom, width, height]
        # ax = fig.add_axes(rect_scatter)
        # ax.imshow(frame_rgb)
        # ax.axis("off")
        # #

        points_in_disks = monte_carlo(
            num_samples, frame_rgb, frame_width, frame_height, frame_max
        )

        # # here
        # for _ in range(10000):
        #     x = np.random.uniform(0, frame_width - (frame_width - frame_max))
        #     y = np.random.uniform(0, frame_height)
        #     color = "k" if in_disk(x, y, frame_rgb) else "r"
        #     ax.scatter(x, y, c=color, s=1)
        # plt.show()
        # #

        void_ratio_arr[i] = (num_samples - points_in_disks) / num_samples
        if i % 10 == 0:
            print("Iteration {} / {}".format(i, num_iteration))

    return void_ratio_arr


def compute_void_ratio(frame, frame_width, frame_height, num_samples):
    # Convert the frame to RGB
    frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

    # find where the plate is located to limit the monte carlo sim
    frame_max = find_first_uniform_column(frame_rgb)
    if frame_max == -1:
        print("It is not picking up the proper column")
        print(
            "This is a warning that the code is properly not doing the correct thing."
        )

    # Monte Carlo simulation to estimate the ratio of disk to total area
    void_ratio_arr = outer_loop(
        num_samples, frame_rgb, frame_width, frame_height, frame_max
    )
    void_ratio_mean = np.mean(void_ratio_arr)
    void_ratio_std = np.std(void_ratio_arr)

    return void_ratio_mean, void_ratio_std


def process_video(cap, num_samples):
    # Process each frame
    frame_idx = 0
    vr_means = []
    vr_stds = []

    frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break

        vr_mean, vr_std = compute_void_ratio(
            frame, frame_width, frame_height, num_samples
        )
        vr_means.append(vr_mean)
        vr_stds.append(vr_std)

        print(f"Frame {frame_idx}: Estimated ratio = {vr_mean:.4f}")

        frame_idx += 1

    cap.release()
    return vr_means, vr_stds

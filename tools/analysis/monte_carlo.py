# --------------------------------------
# import modules
# --------------------------------------

import numpy as np
import cv2
from numba import njit, prange

# ------------------------------------------------
# define the functions
# ------------------------------------------------


# Function to check if a point is in a disk
@njit
def in_disk(x, y, image):
    x, y = int(x), int(y)
    color = image[y, x]
    lower_blue = np.array([0, 0, 0])
    upper_blue = np.array([255, 255, 255])
    return np.all(lower_blue <= color) and np.all(color < upper_blue)


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
    points_in_disks = 0
    for _ in prange(num_samples):
        x = np.random.uniform(0, frame_width - (frame_width - frame_max))
        y = np.random.uniform(0, frame_height)
        if in_disk(x, y, frame_rgb):
            points_in_disks += 1
    return points_in_disks


def compute_void_ratio(frame, frame_width, frame_height, num_samples):
    # Convert the frame to RGB
    frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

    # find where the plate is located to limit the monte carlo sim
    frame_max = find_first_uniform_column(frame_rgb)
    if frame_max == -1:
        print("It is not picking up the proper column")
        exit()

    # Monte Carlo simulation to estimate the ratio of disk to total area
    points_in_disks = monte_carlo(num_samples, frame_rgb, frame_width, frame_height, frame_max)
    void_ratio = (num_samples - points_in_disks) / num_samples
    return void_ratio

def process_video(cap, num_samples):
    # Process each frame
    frame_idx = 0
    void_ratios = []

    frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    while cap.isOpened():
        ret, frame = cap.read()
        if not ret:
            break

        void_ratio = compute_void_ratio(frame, frame_width, frame_height, num_samples)
        void_ratios.append(void_ratio)

        if frame_idx % 50 == 0:
            print(f"Frame {frame_idx}: Estimated ratio = {void_ratio:.4f}")

        frame_idx += 1

    cap.release()
    return void_ratios

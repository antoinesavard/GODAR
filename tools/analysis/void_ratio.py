import numpy as np
import matplotlib.pyplot as plt
import tools.utils.files as ff
import os
import sys
import cv2
from numba import njit, prange, set_num_threads, get_num_threads

# Set the number of threads
#set_num_threads(4)
print(f"Using {get_num_threads()} threads for parallel execution\n")

# loading the relevant data
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
try:
    video_path = str(sys.argv[3])
    print("Path of video = {}".format(video_path))
except:
    print("Nor argument provided by sys.")
    video_path = str(input("video_path = "))

# reading the files
print("Reading the files...")    
#filesx = ff.list_files(output_dir, "x", expno)
#filesy = ff.list_files(output_dir, "y", expno)
#filesr = ff.list_files(output_dir, "r", expno)

#x, y, r = (
#    ff.multiload(output_dir, filesx),
#    ff.multiload(output_dir, filesy),
#    ff.multiload(output_dir, filesr),
#)

#x = x[::compression] / sf
#y = y[::compression] / sf
#r = r[::compression] / sf
    
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
def monte_carlo(num_samples, frame_width, framheight, frame_rgb, frame_max):
    points_in_disks = 0
    for _ in prange(num_samples):
        x = np.random.uniform(0, frame_width - (frame_width - frame_max))
        y = np.random.uniform(0, frame_height)
        if in_disk(x, y, frame_rgb):
            points_in_disks += 1
    return points_in_disks

# Number of Monte Carlo samples
num_samples = 10000000

# Open the video file
cap = cv2.VideoCapture(video_path)

if not cap.isOpened():
    print("Error: Could not open video.")
    exit()

frame_count = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

# Process each frame
frame_idx = 0
void_ratios = []

while cap.isOpened():
    ret, frame = cap.read()
    if not ret:
        break

    # Convert the frame to RGB
    frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

    # find where the plate is located to limit the spread in the monte carlo sim
    frame_max = find_first_uniform_column(frame_rgb)
    if frame_max == -1:
        print("It is not picking up the proper column")
        exit()

    # Monte Carlo simulation to estimate the ratio of blue disk area to total area
    points_in_disks = monte_carlo(num_samples, frame_width, frame_height, frame_rgb, frame_max)
    void_ratio = (num_samples - points_in_disks) / num_samples
    void_ratios.append(void_ratio)

    print(f"Frame {frame_idx}: Estimated ratio = {void_ratio:.4f}")

    # Optional: Visualize the frame with sampled points (only for a subset of points)
    #if frame_idx % 100 == 0:  # Visualize every 10th frame
    #    plt.imshow(frame_rgb)
    #    plt.axis('off')
    #    plt.title(f"Frame {frame_idx + 1} with {num_samples} Monte Carlo Samples")
    #    for _ in range(10000):  # Plot only a subset of points for visibility
    #        x = np.random.uniform(0, frame_width - (frame_width - frame_max))
    #        y = np.random.uniform(0, frame_height)
    #        color = 'k' if in_disk(x, y, frame_rgb) else 'r'
    #        plt.scatter(x, y, c=color, s=1)
    #    plt.show()

    frame_idx += 1

cap.release()

# Plot the estimated ratios over time
plt.plot(void_ratios)
plt.title("Estimated void ratio of the between ice floes over time")
plt.xlabel("Frame number")
plt.ylabel("Void ratio")
plt.show()

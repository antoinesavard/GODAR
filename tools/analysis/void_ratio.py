import numpy as np
import tools.utils.files as ff
import os
import sys
import cv2


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
    print("video_path = {}".format(video_path))
except:
    print("No argument provided by sys.")
    video_path = str(input("video_path ="))
    
print("Reading the files...")

filesx = ff.list_files(output_dir, "x", expno)
filesy = ff.list_files(output_dir, "y", expno)
filesr = ff.list_files(output_dir, "r", expno)

x, y, r = (
    ff.multiload(output_dir, filesx),
    ff.multiload(output_dir, filesy),
    ff.multiload(output_dir, filesr),
)

x = x[::compression] / sf
y = y[::compression] / sf
r = r[::compression] / sf

# Function to check if a point is in a blue disk
def in_disk(x, y, image):
    x, y = int(x), int(y)
    color = image[y, x]
    lower_blue = np.array([0, 0, 0])
    upper_blue = np.array([250, 250, 250])
    return np.all(lower_blue <= color) and np.all(color <= upper_blue)

# Number of Monte Carlo samples
num_samples = 100000

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
ratios = []

while cap.isOpened():
    ret, frame = cap.read()
    if not ret:
        break

    # Convert the frame to RGB
    frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

    # Monte Carlo simulation to estimate the ratio of blue disk area to total area
    points_in_disks = 0

    for _ in range(num_samples):
        x = np.random.uniform(0, frame_width)
        y = np.random.uniform(0, frame_height)
        if is_point_in_blue_disk(x, y, frame_rgb):
            points_in_blue += 1

    ratio = points_in_blue / num_samples
    ratios.append(ratio)

    print(f"Frame {frame_idx}: Estimated ratio = {ratio:.4f}")

    # Optional: Visualize the frame with sampled points (only for a subset of points)
    if frame_idx % 10 == 0:  # Visualize every 10th frame
        plt.imshow(frame_rgb)
        plt.title(f"Frame {frame_idx} with Monte Carlo Samples")
        for _ in range(1000):  # Plot only a subset of points for visibility
            x = np.random.uniform(0, frame_width)
            y = np.random.uniform(0, frame_height)
            color = 'b' if is_point_in_blue_disk(x, y, frame_rgb) else 'r'
            plt.scatter(x, y, c=color, s=1)
        plt.show()

    frame_idx += 1

cap.release()

# Plot the estimated ratios over time
plt.plot(ratios)
plt.title("Estimated Ratio of Blue Disk Area to Total Area Over Time")
plt.xlabel("Frame")
plt.ylabel("Ratio")
plt.show()

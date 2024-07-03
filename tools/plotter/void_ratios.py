import numpy as np
import matplotlib.pyplot as plt
import sys
import cv2
import tools.analysis.monte_carlo as tam
import tools.plotter.figures as tpf
import tools.utils.files as tuf
from numba import set_num_threads, get_num_threads


# --------------------------------------
# loading the relevant data
# --------------------------------------

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
    print("No argument provided by sys.")
    video_path = str(input("video_path = "))
try:
    num_threads = int(sys.argv[4])
    print("Number of threads to use: {}".format(num_threads))
except:
    print("No argument provided by sys.")
    num_threads = int(input("num_thread = "))
try:
    num_samples = int(sys.argv[5])
    print("Number of Monte Carlo samples: {}".format(num_samples))
except:
    print("No argument provided by sys.")
    num_samples = int(input("num_samples = "))

    
# Set the number of threads
set_num_threads(num_threads)
print(f"Using {get_num_threads()} threads for parallel execution\n")

# reading the files
# print("Reading the files...")
# filesx = tuf.list_files(output_dir, "x", expno)
# filesy = tuf.list_files(output_dir, "y", expno)
# filesr = tuf.list_files(output_dir, "r", expno)

# x, y, r = (
#    tuf.multiload(output_dir, filesx),
#    tuf.multiload(output_dir, filesy),
#    tuf.multiload(output_dir, filesr),
# )

# x = x[::compression] / sf
# y = y[::compression] / sf
# r = r[::compression] / sf

# ------------------------------------------------
# define some global variables
# ------------------------------------------------

# Open the video file
cap = cv2.VideoCapture(video_path)

if not cap.isOpened():
    print("Error: Could not open video.")
    exit()

void_ratios_means, void_ratios_stds = tam.process_video(cap, num_samples)

#-----------------------------------------------------
# plot the void ratios
#-----------------------------------------------------

tpf.void_ratio_plot(void_ratios_means, void_ratios_stds, expno)

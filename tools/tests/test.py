import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.spatial import Voronoi, voronoi_plot_2d

sys.path.append("/Users/asavard/Documents/GODAR/")
import tools.plotter.figures as tpf
import tools.utils.files as tuf

# ----------------------------------------------------------------------

xaxis_limits = 30  # in km
yaxis_limits = 30  # in km
sf = 1e3  # conversion ratio m <-> km
compression = 1  # data compression

# ----------------------------------------------------------------------

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

print("Reading the files...")

# listing the files to read
filesx = tuf.list_files(output_dir, "x", expno)
filesy = tuf.list_files(output_dir, "y", expno)
filesr = tuf.list_files(output_dir, "r", expno)
filestsigxx = tuf.list_files(output_dir, "tsigxx", expno)
filestsigyy = tuf.list_files(output_dir, "tsigyy", expno)
filestsigxy = tuf.list_files(output_dir, "tsigxy", expno)
filestsigyx = tuf.list_files(output_dir, "tsigyx", expno)

# loading the files in memory
x, y, r, tsigxx, tsigyy, tsigxy, tsigyx = (
    tuf.multiload(output_dir, filesx),
    tuf.multiload(output_dir, filesy),
    tuf.multiload(output_dir, filesr),
    tuf.multiload(output_dir, filestsigxx),
    tuf.multiload(output_dir, filestsigyy),
    tuf.multiload(output_dir, filestsigxy),
    tuf.multiload(output_dir, filestsigyx),
)

# compressing the files
x = x[::compression] / sf
y = y[::compression] / sf
r = r[::compression] / sf
tsigxx = tsigxx[::compression] / sf
tsigyy = tsigyy[::compression] / sf
tsigxy = tsigxy[::compression] / sf
tsigyx = tsigyx[::compression] / sf


# convert it in real stress
tsigxx = tsigxx / r**2 / np.pi
tsigyy = tsigyy / r**2 / np.pi
tsigxy = tsigxy / r**2 / np.pi
tsigyx = tsigyx / r**2 / np.pi

# check dimensions of the data
x = tuf.check_dim(x)
y = tuf.check_dim(y)
r = tuf.check_dim(r)
tsigxx = tuf.check_dim(tsigxx)
tsigyy = tuf.check_dim(tsigyy)
tsigxy = tuf.check_dim(tsigxy)
tsigyx = tuf.check_dim(tsigyx)


# Generate Voronoi plot
def plot_stress_field(x, y, stress):
    points = np.column_stack((x, y))
    vor = Voronoi(points)

    # Create a normalization object for consistent color mapping
    norm = Normalize(vmin=np.min(stress), vmax=np.max(stress))

    # Plot Voronoi tessellation
    fig, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False)

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = stress[
                point_idx
            ]  # Get the corresponding stress value for this region
            plt.fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Scatter plot of particles with the same normalization
    sc = plt.scatter(x, y, c=stress, cmap="viridis", edgecolor="k", norm=norm)
    plt.colorbar(sc, label="Stress")
    plt.title("Stress Distribution in 2D DEM")
    plt.show()


# Example usage for timestep 0
plot_stress_field(x[-1], y[-1], tsigxx[-1])

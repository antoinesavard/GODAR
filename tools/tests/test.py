import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from scipy.spatial import Voronoi, voronoi_plot_2d

sys.path.append("/storage/asavard/DEM/")
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
filestfx = tuf.list_files(output_dir, "tfx", expno)
filestfy = tuf.list_files(output_dir, "tfy", expno)

# loading the files in memory
x, y, r, tsigxx, tsigyy, tsigxy, tsigyx, tfx, tfy = (
    tuf.multiload(output_dir, filesx),
    tuf.multiload(output_dir, filesy),
    tuf.multiload(output_dir, filesr),
    tuf.multiload(output_dir, filestsigxx),
    tuf.multiload(output_dir, filestsigyy),
    tuf.multiload(output_dir, filestsigxy),
    tuf.multiload(output_dir, filestsigyx),
    tuf.multiload(output_dir, filestfx),
    tuf.multiload(output_dir, filestfy),
)

# compressing the files
x = x[::compression] / sf
y = y[::compression] / sf
r = r[::compression] / sf
tsigxx = tsigxx[::compression] / sf
tsigyy = tsigyy[::compression] / sf
tsigxy = tsigxy[::compression] / sf
tsigyx = tsigyx[::compression] / sf
tfx = tfx[::compression] / sf
tfy = tfy[::compression] / sf

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
tfx = tuf.check_dim(tfx)
tfy = tuf.check_dim(tfy)

sigma = np.stack((tsigxx,tsigyy,tsigxy,tsigyx),axis=-1)
tf = np.stack((tfx, tfy), axis=-1)

# Generate Voronoi plot
def plot_stress_field(x, y, stress, time):
    x = x[time-1]
    y = y[time-1]
    stress = stress[time-1]
    points = np.column_stack((x, y))
    vor = Voronoi(points)

    # Create a normalization object for consistent color mapping
    norm = Normalize(vmin=-5e7, vmax=5e7)
    im = cm.ScalarMappable(norm=norm)
    
    # Plot Voronoi tessellation
    fig, ax = plt.subplots(2,2)
    voronoi_plot_2d(vor, ax=ax[0,0], show_points=False, show_vertices=False)
    voronoi_plot_2d(vor, ax=ax[0,1], show_points=False, show_vertices=False)
    voronoi_plot_2d(vor, ax=ax[1,0], show_points=False, show_vertices=False)
    voronoi_plot_2d(vor, ax=ax[1,1], show_points=False, show_vertices=False)

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = stress[point_idx,0]
            # Get the corresponding stress value for this region
            ax[0,0].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = stress[point_idx,1]
            # Get the corresponding stress value for this region
            ax[1,1].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = stress[point_idx,2]
            # Get the corresponding stress value for this region
            ax[0,1].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))


    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = stress[point_idx,3]
            # Get the corresponding stress value for this region
            ax[1,0].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

            
    # Scatter plot of particles with the same normalization
    ax[0,0].scatter(x, y, c=stress[...,0], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    ax[0,1].scatter(x, y, c=stress[...,2], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    ax[1,0].scatter(x, y, c=stress[...,3], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    ax[1,1].scatter(x, y, c=stress[...,1], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    fig.colorbar(im, label="Nm$^{-2}$", ax=ax.ravel().tolist())
    for axis in ax.flat:
        axis.set_xlim(0, xaxis_limits-13)
        axis.set_ylim(0, yaxis_limits)
    ax[0,0].title.set_text("$\sigma_{xx}$")
    ax[1,1].title.set_text("$\sigma_{yy}$")
    ax[0,1].title.set_text("$\sigma_{xy}$")
    ax[1,0].title.set_text("$\sigma_{yx}$")
    fig.suptitle("$\sigma_{ij}$ at $t=$"+str(time))
    fig.savefig("../plots/plot/stress_voronoi_{}.pdf".format(time))
    

# Generate Voronoi plot
def plot_force_field(x, y, force, time):
    x = x[time-1]
    y = y[time-1]
    force = force[time-1]
    points = np.column_stack((x, y))
    vor = Voronoi(points)

    # Create a normalization object for consistent color mapping
    norm = Normalize(vmin=-5e4, vmax=5e4)
    im = cm.ScalarMappable(norm=norm)
    
    # Plot Voronoi tessellation
    fig, ax = plt.subplots(1,2)
    voronoi_plot_2d(vor, ax=ax[0], show_points=False, show_vertices=False)
    voronoi_plot_2d(vor, ax=ax[1], show_points=False, show_vertices=False)

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = force[point_idx,0]
            # Get the corresponding stress value for this region
            ax[0].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = force[point_idx,1]
            # Get the corresponding stress value for this region
            ax[1].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Scatter plot of particles with the same normalization
    ax[0].scatter(x, y, c=force[...,0], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    ax[1].scatter(x, y, c=force[...,1], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    fig.colorbar(im, label="N", ax=ax.ravel().tolist())
    for axis in ax.flat:
        axis.set_xlim(0, xaxis_limits-13)
        axis.set_ylim(0, yaxis_limits)
    ax[0].title.set_text("$F_x$")
    ax[1].title.set_text("$F_y$")
    fig.suptitle("$F_i$ at $t=$"+str(time))
    fig.savefig("../plots/plot/force_voronoi_{}.pdf".format(time))


def plot_force_field(x, y, force, time):
    x = x[time-1]
    y = y[time-1]
    force = force[time-1]
    points = np.column_stack((x, y))
    vor = Voronoi(points)

    # Create a normalization object for consistent color mapping
    norm = Normalize(vmin=-5e4, vmax=5e4)
    im = cm.ScalarMappable(norm=norm)
    
    # Plot Voronoi tessellation
    fig, ax = plt.subplots(1,2)
    voronoi_plot_2d(vor, ax=ax[0], show_points=False, show_vertices=False)
    voronoi_plot_2d(vor, ax=ax[1], show_points=False, show_vertices=False)

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = force[point_idx,0]
            # Get the corresponding stress value for this region
            ax[0].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = force[point_idx,1]
            # Get the corresponding stress value for this region
            ax[1].fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Scatter plot of particles with the same normalization
    ax[0].scatter(x, y, c=force[...,0], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    ax[1].scatter(x, y, c=force[...,1], cmap="viridis", edgecolor="k", norm=norm, alpha=0)
    fig.colorbar(im, label="N", ax=ax.ravel().tolist())
    for axis in ax.flat:
        axis.set_xlim(0, xaxis_limits-13)
        axis.set_ylim(0, yaxis_limits)
    ax[0].title.set_text("$F_x$")
    ax[1].title.set_text("$F_y$")
    fig.suptitle("$F_i$ at $t=$"+str(time))
    fig.savefig("../plots/plot/force_voronoi_{}.pdf".format(time))

def plot_force_quiver(x, y, force, time):
    x = x[time-1]
    y = y[time-1]
    force = force[time-1]
    force_norm = np.sqrt(force[:,0] ** 2 + force[:,1] ** 2)
    points = np.column_stack((x, y))
    vor = Voronoi(points)

    # Create a normalization object for consistent color mapping
    norm = Normalize(vmin=0, vmax=5e4)
    im = cm.ScalarMappable(norm=norm)
    
    # Plot Voronoi tessellation
    fig, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False)

    # Color the regions based on normalized stress
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not -1 in region and region != []:
            polygon = [vor.vertices[i] for i in region]
            value = force_norm[point_idx]
            # Get the corresponding stress value for this region
            ax.fill(*zip(*polygon), color=plt.cm.viridis(norm(value)))

    # Scatter plot of particles with the same normalization
    ax.quiver(x, y, force[:,0], force[:,1])
    fig.colorbar(im, label="N", ax=ax)
    ax.set_xlim(0, xaxis_limits-13)
    ax.set_ylim(0, yaxis_limits)
    fig.suptitle("$F$ at $t=$"+str(time))
    fig.savefig("../plots/plot/force_quiver_{}.pdf".format(time))

    
plot_stress_field(x, y, sigma, 50)
plot_stress_field(x, y, sigma, 100)

plot_force_field(x, y, tf, 50)
plot_force_field(x, y ,tf, 100)

plot_force_quiver(x, y, tf, 50)
plot_force_quiver(x, y, tf, 100)

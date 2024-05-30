# file to create experience 06 initial conditions with 2500 particles

import numpy as np
import sys

# --------------------------------------------
# all inputs for the run
# --------------------------------------------


def reading_input():
    # approximate number of particles in sim
    try:
        n = int(sys.argv[1])
        print("number of particles = {}".format(n))
    except:
        print("No argument provided by sys.")
        n = int(input("number of particles = "))

    # parameters for the radius
    try:
        param1_r = int(sys.argv[2])
        print("mean radius = {}".format(param1_r))
    except:
        print("No argument provided by sys.")
        param1_r = int(input("mean radius = "))

    try:
        param2_r = int(sys.argv[3])
        print("sigma of radius = {}".format(param2_r))
    except:
        print("No argument provided by sys.")
        param2_r = int(input("sigma of radius = "))

    # parameters for the thickness
    try:
        param1_t = int(sys.argv[4])
        print("mean thickness = {}".format(param1_t))
    except:
        print("No argument provided by sys.")
        param1_t = int(input("mean tickness = "))

    try:
        param2_t = int(sys.argv[5])
        print("sigma of thickness = {}".format(param2_t))
    except:
        print("No argument provided by sys.")
        param2_t = int(input("sigma of thickness = "))

    try:
        offset = int(sys.argv[6])
        print("offset of thickness = {}".format(offset))
    except:
        print("No argument provided by sys.")
        offset = int(input("offset of thickness = "))

    try:
        cutoff = int(sys.argv[7])
        print("cutoff of thickness = {}".format(cutoff))
    except:
        print("No argument provided by sys.")
        cutoff = int(input("cutoff of thickness = "))

    # choice of distribution (radius_thickness)
    try:
        dist = int(sys.argv[8])
        print(
            "distribution to use (radius_thickness) in [log, gauss, uni]: {}".format(
                dist
            )
        )
    except:
        print("No argument provided by sys.")
        dist = int(
            input("distribution to use (radius_thickness) in [log, gauss, uni]: ")
        )

    # unique identifier
    try:
        adn = int(sys.argv[9])
        print("unique identifier: {}".format(adn))
    except:
        print("No argument provided by sys.")
        adn = int(input("unique identifier: "))

    return n, param1_r, param2_r, param1_t, param2_t, offset, cutoff, dist


def init(n):
    # size factor
    sf = 1e3
    disks_num_y = np.sqrt(n).astype(int)
    disks_num_x = disks_num_y.astype(int)
    disks_num = disks_num_y * disks_num_x
    return sf, disks_num_x, disks_num_y, disks_num


# --------------------------------------------
# all possible functions to choose from
# --------------------------------------------


def log_log(mean_r, sigma_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.lognormal(mean_r, sigma_r, disks_num)
    thick = np.random.lognormal(mean_t, sigma_t, disks_num) - offset
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def log_gauss(mean_r, sigma_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.lognormal(mean_r, sigma_r, disks_num)
    thick = np.random.normal(mean_t, sigma_t, disks_num)
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def log_uni(mean_r, sigma_r, low_t, high_t, disks_num, offset, cutoff):
    radius = np.random.lognormal(mean_r, sigma_r, disks_num)
    thick = np.random.uniform(low_t, high_t, disks_num)
    return radius, thick


def gauss_log(mean_r, sigma_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.normal(mean_r, sigma_r, disks_num)
    thick = np.random.lognormal(mean_t, sigma_t, disks_num) - offset
    radius = np.where(radius <= 50, 50, radius)
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def gauss_gauss(mean_r, sigma_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.normal(mean_r, sigma_r, disks_num)
    thick = np.random.normal(mean_t, sigma_t, disks_num)
    radius = np.where(radius <= offset, offset, radius)
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def gauss_uni(mean_r, sigma_r, low_t, high_t, disks_num, offset, cutoff):
    radius = np.random.normal(mean_r, sigma_r, disks_num)
    thick = np.random.uniform(low_t, high_t, disks_num)
    radius = np.where(radius <= offset, offset, radius)
    return radius, thick


def uni_log(low_r, high_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.uniform(low_r, high_r, disks_num)
    thick = np.random.lognormal(mean_t, sigma_t, disks_num) - offset
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def uni_gauss(low_r, high_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.uniform(low_r, high_r, disks_num)
    thick = np.random.normal(mean_t, sigma_t, disks_num)
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def uni_uni(low_r, high_r, low_t, high_t, disks_num, offset, cutoff):
    radius = np.random.uniform(low_r, high_r, disks_num)
    thick = np.random.uniform(low_t, high_t, disks_num)
    return radius, thick


# dictionnary of functions to choose from
functions = dict(
    log_log=log_log,
    log_gauss=log_gauss,
    log_uni=log_uni,
    gauss_log=gauss_log,
    gauss_gauss=gauss_gauss,
    gauss_uni=gauss_uni,
    uni_log=uni_log,
    uni_gauss=uni_gauss,
    uni_uni=uni_uni,
)

# --------------------------------------------
# the creation of the input files
# --------------------------------------------


def file_creation(
    sf,
    disks_num_x,
    disks_num_y,
    disks_num,
    param1_r,
    param2_r,
    param1_t,
    param2_t,
    offset,
    cutoff,
    dist,
    adn,
):
    lines = sf * np.arange(1, 2 * disks_num_y + 1, 2)
    lines = lines.astype(str)
    lines = lines.tolist()
    with open("files/y" + adn + ".dat", "w") as f:
        for line in lines:
            i = 0
            while i < disks_num_x:
                f.write(line)
                f.write("\n")
                i += 1

    lines = sf * np.arange(1, 2 * disks_num_x + 1, 2)
    lines = lines.astype(str)
    lines = lines.tolist()
    with open("files/x" + adn + ".dat", "w") as f:
        i = 0
        while i < disks_num_y:
            for line in lines:
                f.write(line)
                f.write("\n")
            i += 1

    f = functions[dist]
    radius, thick = f(param1_r, param2_r, param1_t, param2_t, disks_num, offset, cutoff)

    lines = radius.astype(str).tolist()
    with open("files/r" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")

    lines = thick.astype(str).tolist()
    with open("files/h" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")

    other = np.zeros_like(radius)
    lines = other.astype(str).tolist()
    with open("files/theta" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")

    with open("files/omega" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")


# -----------------------------------------
# bulk of the code
# -----------------------------------------

try:
    print("Reading the init file.")
    with open("plots/input_init.dat", "r") as f:
        lines = f.read().split()
        n = np.asarray(lines[0::9]).astype(int)
        param1_r = np.asarray(lines[1::9]).astype(float)
        param2_r = np.asarray(lines[2::9]).astype(float)
        param1_t = np.asarray(lines[3::9]).astype(float)
        param2_t = np.asarray(lines[4::9]).astype(float)
        offset = np.asarray(lines[5::9]).astype(float)
        cutoff = np.asarray(lines[6::9]).astype(float)
        dist = np.asarray(lines[7::9]).astype(str)
        adn = np.asarray(lines[8::9]).astype(str)
    print("Data written.")

except:
    print("The input file for initiation of particles does not exist.")
    print("Or maybe you are not running this script from the base directory /GODAR")
    n, param1_r, param2_r, param1_t, param2_t, offset, cutoff, dist = reading_input()

for i in range(len(n)):
    sf, disks_num_x, disks_num_y, disks_num = init(n[i])
    file_creation(
        sf,
        disks_num_x,
        disks_num_y,
        disks_num,
        param1_r[i],
        param2_r[i],
        param1_t[i],
        param2_t[i],
        offset[i],
        cutoff[i],
        dist[i],
        adn[i],
    )

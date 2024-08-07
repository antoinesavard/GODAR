# file to create experience 06 initial conditions with 2500 particles

import numpy as np
import sys
import os

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

    # number of particles required in y with the appropriate spacing
    try:
        disks_num_y = int(sys.argv[2])
        print("number of particles in y = {}".format(disks_num_y))
    except:
        print("No argument provided by sys.")
        disks_num_y = int(input("number of particles = "))

    # parameters for the radius
    try:
        param1_r = int(sys.argv[3])
        print("mean radius = {}".format(param1_r))
    except:
        print("No argument provided by sys.")
        param1_r = int(input("mean radius = "))

    try:
        param2_r = int(sys.argv[4])
        print("sigma of radius = {}".format(param2_r))
    except:
        print("No argument provided by sys.")
        param2_r = int(input("sigma of radius = "))

    # parameters for the thickness
    try:
        param1_t = int(sys.argv[5])
        print("mean thickness = {}".format(param1_t))
    except:
        print("No argument provided by sys.")
        param1_t = int(input("mean tickness = "))

    try:
        param2_t = int(sys.argv[6])
        print("sigma of thickness = {}".format(param2_t))
    except:
        print("No argument provided by sys.")
        param2_t = int(input("sigma of thickness = "))

    try:
        offset = int(sys.argv[7])
        print("offset of thickness = {}".format(offset))
    except:
        print("No argument provided by sys.")
        offset = int(input("offset of thickness = "))

    try:
        cutoff = int(sys.argv[8])
        print("cutoff of thickness = {}".format(cutoff))
    except:
        print("No argument provided by sys.")
        cutoff = int(input("cutoff of thickness = "))

    # choice of distribution (radius_thickness)
    try:
        dist = int(sys.argv[9])
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
        adn = int(sys.argv[10])
        print("unique identifier: {}".format(adn))
    except:
        print("No argument provided by sys.")
        adn = int(input("unique identifier: "))

    return n, param1_r, param2_r, param1_t, param2_t, offset, cutoff, dist, adn


def init(n, disks_num_y):
    # size factor
    sf = 1e3
    disks_num_x = (n / disks_num_y).astype(int)
    disks_num = disks_num_y * disks_num_x
    return sf, disks_num_x, disks_num


# --------------------------------------------
# all possible functions to choose from
# --------------------------------------------


def log_log(mean_r, sigma_r, mean_t, sigma_t, disks_num, offset, cutoff):
    radius = np.random.lognormal(mean_r, sigma_r, disks_num)
    thick = np.random.lognormal(mean_t, sigma_t, disks_num) - offset
    thick = np.where(thick <= cutoff, cutoff, thick)
    return radius, thick


def log_const(mean_r, sigma_r, value, void, disks_num, offset, cutoff):
    radius = np.random.lognormal(mean_r, sigma_r, disks_num)
    thick = np.ones(disks_num) * value
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


def gauss_const(mean_r, sigma_r, value, void, disks_num, offset, cutoff):
    radius = np.random.normal(mean_r, sigma_r, disks_num)
    thick = np.ones(disks_num) * value
    radius = np.where(radius <= offset, offset, radius)
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


def uni_const(low_r, high_r, value, void, disks_num, offset, cutoff):
    radius = np.random.uniform(low_r, high_r, disks_num)
    thick = np.ones(disks_num) * value
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


def const_const(value_r, sigma_r, value_h, void, disks_num, offset, cutoff):
    radius = np.ones(disks_num) * value_r
    thick = np.ones(disks_num) * value_h
    return radius, thick


# dictionnary of functions to choose from
functions = dict(
    log_log=log_log,
    log_const=log_const,
    log_gauss=log_gauss,
    log_uni=log_uni,
    gauss_log=gauss_log,
    gauss_const=gauss_const,
    gauss_gauss=gauss_gauss,
    gauss_uni=gauss_uni,
    uni_log=uni_log,
    uni_const=uni_const,
    uni_gauss=uni_gauss,
    uni_uni=uni_uni,
    const_const=const_const,
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
    work_dir,
):
    wall = sf * np.arange(0.5, 2 * disks_num_y, 1)
    walls = np.concatenate((wall, wall))
    walls = walls.astype(str).tolist()
    lines = sf * np.arange(1, 2 * disks_num_y + 1, 2)
    with open(work_dir + "/../files/y" + adn + ".dat", "w") as f:
        for clean_line in lines:
            i = 0
            while i < disks_num_x:
                line = clean_line + np.random.uniform(-1, 1)
                f.write(str(line))
                f.write("\n")
                i += 1
        for line in walls:
            f.write(line)
            f.write("\n")

    wallf = sf * np.ones(2 * disks_num_y) / 2
    wallb = sf * (np.ones_like(wallf) * 2 * disks_num_x + 3 / 2)
    walls = np.concatenate((wallf, wallb))
    walls = walls.astype(str).tolist()
    lines = sf * np.arange(2, 2 * disks_num_x + 1, 2)
    lines = lines.astype(str).tolist()
    with open(work_dir + "/../files/x" + adn + ".dat", "w") as f:
        i = 0
        while i < disks_num_y:
            for line in lines:
                f.write(line)
                f.write("\n")
            i += 1
        for line in walls:
            f.write(line)
            f.write("\n")

    f = functions[dist]
    radius, thick = f(param1_r, param2_r, param1_t, param2_t, disks_num, offset, cutoff)

    print(
        "You should use {} as a search radius in your namelist".format(
            int(np.ceil(np.max(radius)))
        )
    )
    wall = 500 * np.ones(2 * disks_num_y)
    walls = np.concatenate((wall, wall))
    walls = walls.astype(str).tolist()
    lines = radius.astype(str).tolist()
    with open(work_dir + "/../files/r" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")
        for line in walls:
            f.write(line)
            f.write("\n")

    wall = np.ones(2 * disks_num_y)
    walls = np.concatenate((wall, wall))
    walls = walls.astype(str).tolist()
    lines = thick.astype(str).tolist()
    with open(work_dir + "/../files/h" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")
        for line in walls:
            f.write(line)
            f.write("\n")

    other = np.zeros_like(np.concatenate((radius, wall, wall)))
    lines = other.astype(str).tolist()
    with open(work_dir + "/../files/theta" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")

    with open(work_dir + "/../files/omega" + adn + ".dat", "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")

    return len(other)


# -----------------------------------------
# bulk of the code
# -----------------------------------------

try:
    print("Reading the input file.")
    with open("init_args.dat", "r") as f:
        lines = f.read().split()
        first_line_skip = 10
        n = np.asarray(lines[0 + first_line_skip :: 10]).astype(int)
        disks_num_y = np.asarray(lines[1 + first_line_skip :: 10]).astype(int)
        param1_r = np.asarray(lines[2 + first_line_skip :: 10]).astype(float)
        param2_r = np.asarray(lines[3 + first_line_skip :: 10]).astype(float)
        param1_t = np.asarray(lines[4 + first_line_skip :: 10]).astype(float)
        param2_t = np.asarray(lines[5 + first_line_skip :: 10]).astype(float)
        offset = np.asarray(lines[6 + first_line_skip :: 10]).astype(float)
        cutoff = np.asarray(lines[7 + first_line_skip :: 10]).astype(float)
        dist = np.asarray(lines[8 + first_line_skip :: 10]).astype(str)
        adn = np.asarray(lines[9 + first_line_skip :: 10]).astype(str)
    work_dir = os.getcwd()

except:
    print("The input file for initiation of particles was not provided.")
    print("Or maybe you are running this script from /GODAR/tools/init/")
    print("In that case, I will read inputs that you provide me with.")

    # we probably are in /GODAR/tools/init/
    path = os.getcwd()
    os.chdir(path + "/..")
    work_dir = os.getcwd()

    (
        n,
        disks_num_y,
        param1_r,
        param2_r,
        param1_t,
        param2_t,
        offset,
        cutoff,
        dist,
        adn,
    ) = reading_input()

for i in range(len(n)):
    sf, disks_num_x, disks_num = init(n[i], disks_num_y[i])
    total_num = file_creation(
        sf,
        disks_num_x,
        disks_num_y[i],
        disks_num,
        param1_r[i],
        param2_r[i],
        param1_t[i],
        param2_t[i],
        offset[i],
        cutoff[i],
        dist[i],
        adn[i],
        work_dir,
    )

print(
    "You will need {} namelists when running the duplication.sh program".format(len(n))
)
print("Compile godar with {} particles".format(total_num))
print("Compile godar with {} km in x".format((disks_num_x + 1) * 2))
print("Compile godar with {} km in y".format(disks_num_y * 2))

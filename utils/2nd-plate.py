# file to create experience 07 initial conditions with 5000 particles

import numpy as np
import sys

gap = 30000


def reading_input():
    # unique identifier
    try:
        adn = int(sys.argv[1])
        print("unique identifier: {}".format(adn))
    except:
        print("No argument provided by sys.")
        adn = int(input("unique identifier: "))

    # approximate number of particles in sim
    try:
        n = int(sys.argv[2])
        print("number of particles = {}".format(n))
    except:
        print("No argument provided by sys.")
        n = int(input("number of particles = "))

    return adn, n


try:
    print("Reading the pyargs.dat file.")
    with open("pyargs.dat", "r") as f:
        lines = f.read().split()
        adn = np.asarray(lines[0::2]).astype(int)
        n = np.asarray(lines[1::2]).astype(int)

except:
    print("The input file for initiation of particles does not exist.")
    print("Or maybe you are not running this script from /GODAR/utils")
    adn, n = reading_input()

adn_new = adn + len(adn)

for j in range(len(adn)):
    # x
    with open("../output/x." + str(adn[j]), "r") as f:
        xlines1 = f.read().split()
        xlines1 = np.asarray(xlines1).astype(float).reshape(1000, n[j])
        xlines2 = xlines1[-1] + gap
        xlines = np.concatenate((xlines1[-1], xlines2))

    with open("../files/x" + str(adn_new[j]) + ".dat", "w") as f:
        for i, xline in enumerate(xlines):
            f.write(str(xline))
            f.write("\n")

    # y
    with open("../output/y." + str(adn[j]), "r") as f:
        ylines1 = f.read().split()
        ylines1 = np.asarray(ylines1).astype(float).reshape(1000, n[j])
        ylines2 = ylines1[-1]
        ylines = np.concatenate((ylines1[-1], ylines2))

    with open("../files/y" + str(adn_new[j]) + ".dat", "w") as f:
        for i, yline in enumerate(ylines):
            f.write(str(yline))
            f.write("\n")

    # r
    with open("../output/r." + str(adn[j]), "r") as f:
        rlines1 = f.read().split()
        rlines1 = np.asarray(rlines1).astype(float).reshape(1000, n[j])
        rlines2 = rlines1[-1]
        rlines = np.concatenate((rlines1[-1], rlines2))

    with open("../files/r" + str(adn_new[j]) + ".dat", "w") as f:
        for i, rline in enumerate(rlines):
            f.write(str(rline))
            f.write("\n")

    # h
    with open("../output/h." + str(adn[j]), "r") as f:
        hlines1 = f.read().split()
        hlines1 = np.asarray(hlines1).astype(float).reshape(1000, n[j])
        hlines2 = hlines1[-1]
        hlines = np.concatenate((hlines1[-1], hlines2))

    with open("../files/h" + str(adn_new[j]) + ".dat", "w") as f:
        for i, hline in enumerate(hlines):
            f.write(str(hline))
            f.write("\n")

    # theta
    with open("../output/theta." + str(adn[j]), "r") as f:
        tlines1 = f.read().split()
        tlines1 = np.asarray(tlines1).astype(float).reshape(1000, n[j])
        tlines2 = tlines1[-1]
        tlines = np.concatenate((tlines1[-1], tlines2))

    with open("../files/theta" + str(adn_new[j]) + ".dat", "w") as f:
        for i, tline in enumerate(tlines):
            f.write(str(tline))
            f.write("\n")

    # omega
    with open("../output/omega." + str(adn[j]), "r") as f:
        olines1 = f.read().split()
        olines1 = np.asarray(olines1).astype(float).reshape(1000, n[j])
        olines2 = olines1[-1]
        olines = np.concatenate((olines1[-1], olines2))

    with open("../files/omega" + str(adn_new[j]) + ".dat", "w") as f:
        for i, oline in enumerate(olines):
            f.write(str(oline))
            f.write("\n")

    print("Done with line {}/{}".format(j+1,len(adn)))

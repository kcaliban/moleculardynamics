import matplotlib.pyplot as plt
import numpy as np
import sys
import pathlib

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_e_t.py <file>")
        sys.exit(1)

    file = sys.argv[1]
    with open(file, "r") as f:
        lines = f.readlines()

    y = []
    for line in lines:
        if line.startswith("#"):
            continue
        line = line.split(' ')
        y.append(float(line[1]))
        #x.append(float(line[0]))

    x = np.arange(0, len(y), 1)

    poly = np.polyfit(x, y, 10)
    poly_y = np.poly1d(poly)(x)
    #plt.xscale("log")
    #plt.yscale("log")
    plt.plot(x, y, marker="")
    # plt.plot(x, poly_y, marker="")
    #plt.xticks(np.arange(min(x), max(x)+1, 100))
    #plt.xticks(np.arange(-15000, -2000, 1000))
    plt.yticks(np.arange(100, 2000, 100))
    plt.xlabel("total energy (eV)")
    plt.ylabel("temperature (K)")
    plt.savefig("E_t_" + pathlib.Path(file).stem + ".png")

import matplotlib.pyplot as plt
import sys
import pathlib

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_total_energy.py <file>")
        sys.exit(1)

    file = sys.argv[1]
    with open(file, "r") as f:
        lines = f.readlines()

    x = []
    y = []
    for line in lines:
        if line.startswith("#"):
            continue
        line = line.split('\t')
        x.append(float(line[0]))
        y.append(float(line[1]))

    plt.plot(x, y)
    plt.xlabel("time (fs)")
    plt.ylabel("total energy (kJ/mol)")
    plt.savefig("total_energy" + pathlib.Path(file).stem + ".png")

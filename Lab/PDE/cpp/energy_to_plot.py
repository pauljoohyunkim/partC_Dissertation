import json
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: energy_to_plot.py [energy.json]")

    with open(sys.argv[1], "r") as file:
        data = json.load(file)

    plt.plot(data)


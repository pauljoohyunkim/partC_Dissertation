#!/bin/python

import matplotlib.pyplot as plt
import json
import sys

if __name__ == "__main__":
    try:
        with open(sys.argv[1], "r") as file:
            x, us = json.load(file)
    except IndexError:
        print("Usage: json_plot.py [json file of plot data]")
        sys.exit(1)


    # Load the solution at each time step and plot
    for u in us:
        plt.plot(x, u)

    plt.show()

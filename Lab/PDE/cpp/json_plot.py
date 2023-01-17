#!/bin/python

import matplotlib.pyplot as plt
import json
import sys

if __name__ == "__main__":
    # Skip time steps by this number
    skipnum = 1
    try:
        with open(sys.argv[1], "r") as file:
            x, us = json.load(file)
        if len(sys.argv) > 2:
            skipnum = int(sys.argv[2])
    except IndexError:
        print("Usage: json_plot.py [json file of plot data]")
        sys.exit(1)

    try:
        #for u in us:
        #    plt.plot(x, u)
        for i in range(0, len(us), skipnum):
            plt.plot(x, us[i])

        plt.show()
    except KeyboardInterrupt:
        print("Keyboard Interrupt Received. Quitting")
        sys.exit(1)
    
    except:
        print("Error while plotting...")
        sys.exit(1)

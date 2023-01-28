import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: json_to_3D_anim.py [path to x.json] [... to y.json] [... to z.json]")
        sys.exit(1)

    with open(sys.argv[1], "r") as file:
        xdata = np.array(json.load(file))
        file.close()
    with open(sys.argv[2], "r") as file:
        ydata = np.array(json.load(file))
        file.close()
    with open(sys.argv[3], "r") as file:
        zdata = np.array(json.load(file))
        file.close()
    
    M, J = xdata.shape
    
    # Plotting
    fig = plt.figure()
    ax = plt.subplot(projection="3d")
    line, = ax.plot([],[],[], "ro")
    xpoints = np.zeros(J, dtype=float)
    ypoints = np.zeros(J, dtype=float)
    zpoints = np.zeros(J, dtype=float)
    def init():
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)
        ax.set_zlim(-2, 2)
        return line,

    def update(t):
        for i in range(J):
            xpoints[i] = xdata[t][i]
            ypoints[i] = ydata[t][i]
            zpoints[i] = zdata[t][i]
        line.set_data_3d(xpoints, ypoints, zpoints)
        ax.set_title(str(t))
        return line,

    ani = FuncAnimation(fig, update, frames=range(J), init_func=init, blit=True)

    plt.show()
    #update(10)
    
    pass


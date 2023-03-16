import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import sys

FPS=60

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: json_to_3D_anim.py [path to x.json] [... to y.json] [... to z.json]")
        sys.exit(1)

    with open(sys.argv[1], "r") as file:
        xdata = np.array(json.load(file))
        file.close()
    with open(sys.argv[2], "r") as file:
        ydata = np.array(json.load(file))
        file.close()
    
    M, J = xdata.shape
    
    # Plotting
    fig = plt.figure()
    ax = plt.subplot()
    line, = ax.plot([],[],"bo")
    xpoints = np.zeros(J, dtype=float)
    ypoints = np.zeros(J, dtype=float)
    def init():
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        ax.set_title("Evolution")
        return line,

    def update(t):
        for i in range(J):
            xpoints[i] = xdata[t][i]
            ypoints[i] = ydata[t][i]
        line.set_data(xpoints, ypoints)
        progressString = f"Progress: {t} / {M} ({'%.2f'%(t/M * 100)}%)"
        ax.set_title(progressString)
        if(t % 10 == 0):
            print(progressString)
        return line,

    ani = FuncAnimation(fig, update, frames=range(M), init_func=init, blit=False)
    
    if len(sys.argv) == 4:
        print("Writing video as " + sys.argv[3])
        videoWriter = FFMpegWriter(fps=FPS)
        ani.save(sys.argv[3], writer=videoWriter)
        print("Done!")

    plt.show()

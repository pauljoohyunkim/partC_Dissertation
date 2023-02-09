import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import sys

FPS=60

if __name__ == "__main__":
    if len(sys.argv) < 4:
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
    line, = ax.plot([],[],[], "bo")
    xpoints = np.zeros(J, dtype=float)
    ypoints = np.zeros(J, dtype=float)
    zpoints = np.zeros(J, dtype=float)
    def init():
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 0.5)
        ax.set_zlim(-0.5, 0.5)
        ax.set_title("Evolution")
        return line,

    def update(t):
        for i in range(J):
            xpoints[i] = xdata[t][i]
            ypoints[i] = ydata[t][i]
            zpoints[i] = zdata[t][i]
        line.set_data_3d(xpoints, ypoints, zpoints)
        progressString = f"Progress: {t} / {M} ({'%.2f'%(t/M * 100)}%)"
        ax.set_title(progressString)
        if(t % 10 == 0):
            print(progressString)
        return line,

    ani = FuncAnimation(fig, update, frames=range(M), init_func=init, blit=False)
    
    if len(sys.argv) == 5:
        print("Writing video as " + sys.argv[4])
        videoWriter = FFMpegWriter(fps=FPS)
        ani.save(sys.argv[4], writer=videoWriter)
        print("Done!")

    plt.show()


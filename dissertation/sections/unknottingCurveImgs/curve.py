import numpy as np
import matplotlib.pyplot as plt

# Definition of a Curve
class Curve:
    # Constructor
    def __init__(self, list_of_points, closed=False):
        self.list_of_points = list_of_points
        self.closed = closed
        self.J = len(list_of_points)
    
    # Square Bracket Overload
    def __getitem__(self, key):
        return self.list_of_points[key % self.J]
    def __setitem__(self, key, value):
        self.list_of_points[key] = value
    
    # Length
    def __len__(self):
        return self.J
    
    # Curve Addition
    def __add__(self, otherCurve):
        if len(self) == len(otherCurve) and self.closed == otherCurve.closed:
            return Curve([self[i] + otherCurve[i] for i in range(len(self))])

    # Curve Subtraction
    def __sub__(self, otherCurve):
        if len(self) == len(otherCurve) and self.closed == otherCurve.closed:
            return Curve([self[i] - otherCurve[i] for i in range(len(self))])

    # Scalar Multiplication
    def __mul__(self, val):
        return Curve([self[i] * val for i in range(len(self))])
    
    # Curve Length
    def curveLength(self):
        l = 0
        for i in range(self.J - 1):
            edgeLength = np.linalg.norm(self[i + 1] - self[i])
            l += edgeLength ** 2
        if self.closed:
            edgeLength = np.linalg.norm(self[-1] - self[0])
            l += edgeLength ** 2
        return l



def curvePlot(curve, q2d=False, size=(-5,5), title="", ax=plt):
    J = curve.J

    xpoints = []
    ypoints = []
    zpoints = []
    for i in range(J):
        xpoints.append(curve[i][0])
        ypoints.append(curve[i][1])
        zpoints.append(curve[i][2])
    
    if q2d:
        ax.plot(xpoints, ypoints)
    else:
        ax.plot(xpoints, ypoints, zpoints)
        #ax.set_zlim(size)
    #ax.set_xlim(size)
    #ax.set_ylim(size)
    #if title:
        #ax.set_title(title)
    #plt.show()

import numpy as np
from copy import deepcopy

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



# Simple Curve Energy Version 2 (Averaged over 4 points)
def simpleCurveEnergy2(curve):
    J = curve.J

    integral = 0
    for i in range(J):
        for j in range(J):
            # Undefined for |i - j| <= 1
            if abs(i - j) <= 1 or abs(i - j + J) <= 1 or abs(i - j - J) <= 1:
                continue

            contribution = 1 / np.linalg.norm(curve[i] - curve[j]) + 1 / np.linalg.norm(curve[i] - curve[j+1])
            contribution += 1 / np.linalg.norm(curve[i+1] - curve[j]) + 1 / np.linalg.norm(curve[i+1] - curve[j+1])
            contribution /= 4
            integral += contribution * np.linalg.norm(curve[i] - curve[i+1]) * np.linalg.norm(curve[j] - curve[j+1])
    
    return integral

# xa is numpy array of a_n in Fourier Series
# Similarly for others
# It is expected that all arrays passed have the same shape
# resolution is the number of points of the curve.
# Note that the index-0 item of *b lists are ignored
def curveFromFourierCoefficients(xa, xb, ya, yb, za, zb, resolution=10):
    J = len(xa) - 1

    def fourierFromCoefficient(an, bn):
        return lambda x : an[0]/2 + sum([np.cos(x)*an[i] + np.sin(x)*bn[i] for i in range(1, J)])
    
    fx = fourierFromCoefficient(xa, xb)
    fy = fourierFromCoefficient(ya, yb)
    fz = fourierFromCoefficient(za, zb)

    points = []
    for i in range(resolution):
        theta = 2 * np.pi * i / resolution
        points.append(np.array([
            fx(theta),
            fy(theta),
            fz(theta)
        ]))
    
    return Curve(points, True)
import numpy as np
import discrete

# p, q, T are all array
def kernelalphabeta(p, q, T, alpha=2, beta=4):
    pmq = p - q
    numerator = discrete.discretelPNorm(np.cross(T, pmq)) ** alpha
    denominator = discrete.discretelPNorm(pmq) ** beta
    return numerator / denominator

# Pass a list of points.
def energy(points):
    J = len(points)

    e = 0
    for i in range(J):
        for j in range(J):
            if abs(i - j) > 1 and abs(i - j + J) > 1 and abs(i - j - J) > 1:
                xi = points[i]
                xipm = points[i+1]
                xj = points[j]
                xjpm = points[j+1]
                xI = xipm - xi
                lI = discrete.discretelPNorm(xI)
                TI = xI / lI
                e += (kernelalphabeta(xi, xj, TI) + kernelalphabeta(xi, xjpm, TI) + kernelalphabeta(xipm, xj, TI) + kernelalphabeta(xipm, xjpm, TI)) * 0.25 * lI * discrete.discretelPNorm(xjpm-xj)

    return e

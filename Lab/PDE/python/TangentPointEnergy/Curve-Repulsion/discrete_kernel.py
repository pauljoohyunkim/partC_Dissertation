import numpy as np
import discrete

# p, q, T are all array
def kernelalphabeta(p, q, T, alpha=2, beta=4):
    pmq = p - q
    numerator = discrete.discretelPNorm(np.cross(T, pmq), 1) ** alpha
    denominator = discrete.discretelPNorm(pmq, 1) ** beta
    return numerator / denominator

# Pass a list of points.
def energy(points):
    J = len(points)

    e = 0
    for i in range(J):
        for j in range(J):
            if abs(i - j) > 1 and abs(i - j + J) > 1 and abs(i - j - J) > 1:
                xi = points[i % J]
                xipm = points[(i+1) % J]
                xj = points[j % J]
                xjpm = points[(j+1) % J]
                xI = xipm - xi
                lI = discrete.discretelPNorm(xI, 1)
                TI = xI / lI
                kernelAB = 0
                kernelAB += kernelalphabeta(xi, xj, TI) 
                kernelAB += kernelalphabeta(xi, xjpm, TI) 
                kernelAB += kernelalphabeta(xipm, xj, TI) 
                kernelAB += kernelalphabeta(xipm, xjpm, TI)
                kernelAB *= 0.25 * lI * discrete.discretelPNorm(xjpm-xj, 1)
                e += kernelAB
                print(f"{i}, {j}: kernelAB: {kernelAB}")

    return e

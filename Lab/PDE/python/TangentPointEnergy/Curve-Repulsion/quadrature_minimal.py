import numpy as np
from curve import Curve
from tangentPointEnergy import kernelalphabeta
from quadrature4 import dkalphabeta, derivative_index, dProductOfLengths

# k_ij in the sum for energy
def kij(curve, i, j, alpha=2, beta=4):
    TI = curve[i+1] - curve[i]
    TI = TI / np.linalg.norm(TI)
    if np.dot(curve[i+1]-curve[i], curve[j+1]-curve[j]) > 0:
        # "Same Direction"
        if np.linalg.norm(curve[i]-curve[j]) < np.linalg.norm(curve[i+1]-curve[j+1]):
            res = kernelalphabeta(curve[i], curve[j], TI, alpha=alpha, beta=beta)
        else:
            res = kernelalphabeta(curve[i+1], curve[j+1], TI, alpha=alpha, beta=beta)
    else:
        # "Different Direction"
        if np.linalg.norm(curve[i]-curve[j+1]) < np.linalg.norm(curve[i+1]-curve[j]):
            res = kernelalphabeta(curve[i], curve[j+1], TI, alpha=alpha, beta=beta)
        else:
            res = kernelalphabeta(curve[i+1], curve[j], TI, alpha=alpha, beta=beta)
    return res

def dkij(curve, i, j, k, alpha=3, beta=6):
    TI = curve[i+1] - curve[i]
    TI = TI / np.linalg.norm(TI)
    if np.dot(curve[i+1]-curve[i], curve[j+1]-curve[j]) > 0:
        # "Same Direction"
        if np.linalg.norm(curve[i]-curve[j]) < np.linalg.norm(curve[i+1]-curve[j+1]):
            res = dkalphabeta(curve, i, j, i, k, alpha=alpha, beta=beta)
        else:
            res = dkalphabeta(curve, i+1, j+1, i, k, alpha=alpha, beta=beta)
    else:
        # "Different Direction"
        if np.linalg.norm(curve[i]-curve[j+1]) < np.linalg.norm(curve[i+1]-curve[j]):
            res = dkalphabeta(curve, i, j+1, i, k, alpha=alpha, beta=beta)
        else:
            res = dkalphabeta(curve, i+1, j, i, k, alpha=alpha, beta=beta)
    return res

def dEnergy(curve, alpha, beta):
    J = len(curve)
    res = Curve([np.array([0,0,0], dtype=curve[0].dtype) for i in range(J)])

    for k in range(J):
        # Generate indices that are relevant
        dIndex = derivative_index(k, J)
        for i, j in dIndex:
            xiEdgeLen = np.linalg.norm(curve[i+1] - curve[i])
            xjEdgeLen = np.linalg.norm(curve[j+1] - curve[j])
            summand = dkij(curve, i, j, k, alpha, beta) * xiEdgeLen * xjEdgeLen
            summand += kij(curve, i, j, alpha, beta) * dProductOfLengths(curve, i, j, k)
            res[k] += summand
        
    return res


# Tangent-point Curve Energy (One-point Quadrature)
def curveEnergyTangentPointKernelFromMinimal(curve, alpha=3, beta=6):
    J = curve.J

    integral = 0
    for i in range(J):
        for j in range(J):
            # Undefined for |i - j| <= 1
            if abs(i - j) <= 1 or abs(i - j + J) <= 1 or abs(i - j - J) <= 1:
                continue

            # Checking orientation, then picking two shortest points on each edge i and j
            # Check orientation
            xI = curve[i+1] - curve[i]
            contribution = kij(curve, i, j, alpha=alpha, beta=beta)
            integral += contribution * np.linalg.norm(xI) * np.linalg.norm(curve[j] - curve[j+1])
    
    return integral


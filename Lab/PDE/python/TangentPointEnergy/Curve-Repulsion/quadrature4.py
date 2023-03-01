from curve import Curve
from copy import deepcopy
from tangentPointEnergy import kernelalphabeta
import numpy as np


def kij(curve, i, j, alpha=2, beta=4):
    TI = curve[i+1] - curve[i]
    TI = TI / np.linalg.norm(TI)
    res = kernelalphabeta(curve[i], curve[j], TI, alpha, beta)
    res += kernelalphabeta(curve[i], curve[j+1], TI, alpha, beta)
    res += kernelalphabeta(curve[i+1], curve[j], TI, alpha, beta)
    res += kernelalphabeta(curve[i+1], curve[j+1], TI, alpha, beta)
    return res / 4


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
                lI = np.linalg.norm(xI)
                TI = xI / lI
                kernelAB = 0
                kernelAB += kernelalphabeta(xi, xj, TI) 
                kernelAB += kernelalphabeta(xi, xjpm, TI) 
                kernelAB += kernelalphabeta(xipm, xj, TI) 
                kernelAB += kernelalphabeta(xipm, xjpm, TI)
                kernelAB *= 0.25 * lI * np.linalg.norm(xjpm-xj)
                e += kernelAB
                #print(f"{i}, {j}: kernelAB: {kernelAB}")

    return e


# Generates the index pairs that are responsible for the analytic derivative.
# Note that 4(J-1) pairs are generated.
def derivative_index(k, J):
    index_list = []
    for i in range(J):
        if abs(k - i) > 1 and abs(k - i + J) > 1 and abs(k - i - J) > 1:
            index_list.append((k, i))
    for i in range(J):
        if abs(k - 1 - i) > 1 and abs(k - 1 - i + J) > 1 and abs(k - 1 - i - J) > 1:
            index_list.append(((k-1) % J, i))
    for i in range(J):
        if abs(k - i) > 1 and abs(k - i + J) > 1 and abs(k - i - J) > 1:
            index_list.append((i, k))
    for i in range(J):
        if abs(k - 1 - i) > 1 and abs(k - 1 - i + J) > 1 and abs(k - 1 - i - J) > 1:
            index_list.append((i, (k-1) % J))
    

    return index_list

def kjk(curve, p, q, r, alpha=3, beta=6):
    k = p
    j = q
    xkEdge = curve[k+1] - curve[k]
    xkEdgeLen = np.linalg.norm(xkEdge)
    xkj = curve[k] - curve[j]
    xkjLen = np.linalg.norm(xkj)
    xi = xkEdgeLen**2 * xkjLen**2 - (np.dot(xkEdge, xkj))**2
    eta = xkjLen**beta * xkEdgeLen**alpha
    dxi = -2 * xkEdge * xkjLen**2 + 2 * xkEdgeLen**2 * xkj - 2 * np.dot(xkEdge, xkj) * (xkEdge - xkj)
    deta = beta * xkjLen**(beta - 2) * xkEdgeLen**alpha * xkj + alpha * xkjLen**beta * xkEdgeLen**(alpha-2)*(-xkEdge)
    return (xi, eta, dxi, deta)

def ijk(curve, p, q, r, alpha=3, beta=6):
    i = p
    j = q
    k = r
    xkEdge = curve[k+1] - curve[k]
    xkEdgeLen = np.linalg.norm(xkEdge)
    xij = curve[i] - curve[j]
    xijLen = np.linalg.norm(xij)
    xi = xkEdgeLen**2 * xijLen**2 - (np.dot(xkEdge, xij))**2
    eta = xijLen**beta * xkEdgeLen**alpha
    dxi = -2 * xijLen**2 * xkEdge + 2 * np.dot(xkEdge, xij) * xij
    deta = alpha * xijLen**beta * xkEdgeLen**(alpha-2) * (-xkEdge)
    return (xi, eta, dxi, deta)

def km1jkm1(curve, p, q, r, alpha=3, beta=6):
    k = p + 1
    j = q
    xkEdge = curve[k] - curve[k-1]
    xkEdgeLen = np.linalg.norm(xkEdge)
    xkmj = curve[k-1] - curve[j]
    xkmjLen = np.linalg.norm(xkmj)
    xi = xkEdgeLen**2 * xkmjLen**2 - (np.dot(xkEdge, xkmj))**2
    eta = xkmjLen**beta * xkEdgeLen**alpha
    dxi = 2 * xkmjLen**2 * xkEdge - 2 * np.dot(xkEdge, xkmj) * xkmj
    deta = alpha * xkmjLen**beta * xkEdgeLen**(alpha-2) * xkEdge
    return (xi, eta, dxi, deta)

def kjkm1(curve, p, q, r, alpha=3, beta=6):
    k = p
    j = q
    xkEdge = curve[k] - curve[k-1]
    xkEdgeLen = np.linalg.norm(xkEdge)
    xkj = curve[k] - curve[j]
    xkjLen = np.linalg.norm(xkj)
    xi = xkEdgeLen**2 * xkjLen**2 - (np.dot(xkEdge, xkj))**2
    eta = xkjLen**beta * xkEdgeLen**alpha
    dxi = 2 * xkjLen**2 * xkEdge + 2 * xkEdgeLen**2 * xkj - 2 * np.dot(xkEdge, xkj) * (xkEdge + xkj)
    deta = beta * xkEdgeLen**alpha * xkjLen**(beta-2) * xkj + alpha * xkjLen**beta * xkEdgeLen**(alpha-2) * xkEdge
    return (xi, eta, dxi, deta)

def ikj(curve, p, q, r, alpha=3, beta=6):
    i = p
    j = r
    k = q
    xjEdge = curve[j+1] - curve[j]
    xjEdgeLen = np.linalg.norm(xjEdge)
    xki = curve[k] - curve[i]
    xkiLen = np.linalg.norm(xki)
    xi = xjEdgeLen**2 * xkiLen**2 - (np.dot(xjEdge, xki))**2
    eta = xkiLen**beta * xjEdgeLen**alpha
    dxi = 2 * xjEdgeLen**2 * xki - 2 * np.dot(xjEdge, xki) * xjEdge
    deta = beta * xjEdgeLen**alpha * xkiLen**(beta-2) * xki
    return (xi, eta, dxi, deta)

# Derivative of k_{\beta}^{\alpha}
def dkalphabeta(curve, p, q, r, k, alpha=3, beta=6):
    J = len(curve)
    p = p % J
    q = q % J
    r = r % J
    k = k % J
    # (k, j, k)
    if (p, r) == (k, k):
        xiEtaSorter = kjk
        #print("kjk")
    # (i, j, k)
    elif r == k:
        xiEtaSorter = ijk
        #print("ijk")
    elif (p, r) == ((k-1)%J, (k-1)%J):
        xiEtaSorter = km1jkm1
        #print("k-1,j,k-1")
    elif (p, r) == (k, (k-1)%J):
        xiEtaSorter = kjkm1
        #print("k, j, k-1")
    elif q == k:
        xiEtaSorter = ikj
        #print("ikj")
    else:
        #print(f"k={k}: NOT DEFINED: {p}, {q}, {r}")
        return 0.0
    xi, eta, dxi, deta = xiEtaSorter(curve, p, q, r, alpha, beta)

    res = (alpha/2 * xi**(alpha/2 - 1) * dxi * eta - xi**(alpha/2) * deta) / (eta**2)
    return res

def dkij(curve, i, j, k, alpha=3, beta=6):
    res = 0
    res += dkalphabeta(curve, i, j, i, k, alpha, beta)
    res += dkalphabeta(curve, i, j+1, i, k, alpha, beta)
    res += dkalphabeta(curve, i+1, j, i, k, alpha, beta)
    res += dkalphabeta(curve, i+1, j+1, i, k, alpha, beta)
    return res / 4


# Derivative of \norm{x_{i+1} - x_i} \norm{x_{j+1} - x_j}
def dProductOfLengths(curve, p, q, k):
    J = len(curve)
    i = p % J
    j = q % J
    k = k % J
    xiEdge = curve[i+1] - curve[i]
    xiEdgeLen = np.linalg.norm(xiEdge)
    xjEdge = curve[j+1] - curve[j]
    xjEdgeLen = np.linalg.norm(xjEdge)
    if i == k:
        res = xjEdgeLen * (-xiEdge) / xiEdgeLen
    elif (i+1) % J == k:
        res = xjEdgeLen * xiEdge / xiEdgeLen
    elif j == k:
        res = xiEdgeLen * (-xjEdge) / xjEdgeLen
    elif (j+1) % J == k:
        res = xiEdgeLen * xjEdge / xjEdgeLen
    else:
        res = 0
    return res

# Derivative of Energy: Outputs a curve
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




# LEGACY FUNCTIONS
def curveDifferential(curve, index, perturbation_vec, curveEnergy=energy):
    # Perturb the curve in + and - perturb vec.
    # May be improved.
    curvep = Curve(deepcopy(curve.list_of_points))
    curvep[index] += perturbation_vec
    curven = Curve(deepcopy(curve.list_of_points))
    curven[index] -= perturbation_vec

    energyP = curveEnergy(curvep)
    energyN = curveEnergy(curven)
    differential = (energyP - energyN) / 2
    return differential
import numpy as np

# p, q, T are all array
def kernelalphabeta(p, q, T, alpha=2, beta=4):
    pmq = p - q
    numerator = np.linalg.norm(np.cross(T, pmq)) ** alpha
    denominator = np.linalg.norm(pmq) ** beta
    return numerator / denominator
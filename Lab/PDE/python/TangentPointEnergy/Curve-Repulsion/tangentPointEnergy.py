import numpy as np

# p, q, T are all array
def kernelalphabeta(p, q, T, alpha=2, beta=4):
    pmq = p - q
    numerator = np.linalg.norm(np.cross(T, pmq)) ** alpha
    denominator = np.linalg.norm(pmq) ** beta
    return numerator / denominator

def thresholder(x, threshold=1e-10, enabled=False):
    if enabled:
        # Numpy array
        if type(x) == np.ndarray:
            if np.linalg.norm(x, 1) < threshold:
                return np.array([0, 0, 0], dtype=x.dtype)
            else:
                return x
        # Float Value
        else:
            if abs(x) < threshold:
                return 0.0
            else:
                return x

    # Turned off:
    else:
        return x

if __name__ == "__main__":
    f = thresholder(1e-11, enabled=True)

    pass
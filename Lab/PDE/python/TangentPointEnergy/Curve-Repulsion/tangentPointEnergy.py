import numpy as np
import curve

# p, q, T are all array
def kernelalphabeta(p, q, T, alpha=2, beta=4):
    pmq = p - q
    numerator = np.linalg.norm(np.cross(T, pmq)) ** alpha
    denominator = np.linalg.norm(pmq) ** beta
    return numerator / denominator

def thresholder(x, threshold=1e-8, enabled=False):
    if enabled:
        # Numpy array (Threshold each component)
        if isinstance(x, np.ndarray):
            length = len(x)
            res = np.zeros(length, dtype=x.dtype)
            for i in range(length):
                if abs(x[i]) < threshold:
                    res[i] = 0.0
                else:
                    res[i] = x[i]
            return res
        # Float Value
        elif isinstance(x, float):
            if abs(x) < threshold:
                return 0.0
            else:
                return x
        elif isinstance(x, curve.Curve):
            J = len(x)
            for i in range(J):
                x[i] = thresholder(x[i], enabled=True)
            return x

    # Turned off:
    else:
        return x

if __name__ == "__main__":
    f = thresholder(1e-11, enabled=True)

    pass
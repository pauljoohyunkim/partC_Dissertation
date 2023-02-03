

def discretelPNorm(arr, deltaX, p=2):
    norm = 0
    for val in arr:
        norm += abs(val) ** p
    return (norm * deltaX) ** (1/p)

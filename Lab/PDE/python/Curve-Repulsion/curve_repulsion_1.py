import numpy as np
import discrete
import discrete_kernel

if __name__ == "__main__":
    '''
    Points on the curve (or rather... polygon)
    x1(2.5, 3, 2)
    x2(-3.5, 3.5, -2)
    x3(-3.2, -3.25, 1.8)
    x4(5.5, -6.5, 0)
    '''

    '''
    x1 = np.array([2.5, 3, 2])
    x2 = np.array([-3.5, 3.5, -2])
    x3 = np.array([-3.2, -3.25, 1.8])
    x4 = np.array([5.5, -6.5, 0.0])
    '''

    x1 = np.array([1,0,-1], dtype=float)
    x2 = np.array([2,2,-2], dtype=float)
    x3 = np.array([3,4,3], dtype=float)
    x4 = np.array([4,6,-4], dtype=float)
    x5 = np.array([5,8,5], dtype=float)
    x6 = np.array([6,-1,0.6], dtype=float)

    c = [x1, x2, x3, x4, x5, x6]

    discretekernelenergyOrigin = discrete_kernel.energy(c)

    c[5] += np.array([0.0,0,0.1])
    d = [x1, x2, x3, x4, x5, x6]

    discretekernelenergyPerturbed = discrete_kernel.energy(d)

    print(discretekernelenergyOrigin)
    print(discretekernelenergyPerturbed)
    print(discretekernelenergyPerturbed - discretekernelenergyOrigin)


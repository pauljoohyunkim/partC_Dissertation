import numpy as np
import matplotlib.pyplot as plt


# Based on NSPDEs 5.5 Finite Difference Approximation of the Heat Equation
# Parameters
# \{(x,t) \in (a,b) \times [0,T]\}
# J Final spacial point index
# T Final time (in second)
# M Final time step
# deltaT Time mesh size
# deltaX Space mesh size
def heatEvolve(a=-1,b=1,J=20,T=10,M=10000,u_a=0,u_b=0,u_initial = lambda x: 10):
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    cfl = deltaT / (deltaX ** 2)            # Want this below 1/2
    if(cfl > 0.5):
        print(f"Warning: the CFL number is {cfl}, so stability is not guaranteed")
    
    # Explicit Scheme
    u = np.array([[u_a] + [u_initial(a + j * deltaX) for j in range(1,J)] + [u_b]], dtype=float)
    for m in range(M):
        ump1 = np.array([[u_a] + [0 for _ in range(1, J)] + [u_b]], dtype=float)
        # Explicit Scheme
        for j in range(1,J):
            ump1[0,j] = u[m,j] + cfl * (u[m,j+1] - 2 * u[m,j] + u[m,j-1])
        u = np.vstack([u, ump1])


if __name__ == "__main__":
    a = -1
    b = 1
    J = 20              # Final spacial point index
    T = 10              # Final time (in second)
    M = 10000           # Final time step
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    cfl = deltaT / (deltaX ** 2)            # Want this below 1/2
    if(cfl > 0.5):
        print(f"Warning: the CFL number is {cfl}, so stability is not guaranteed")
    # Dirichlet Boundary Conditions u(a) = 0, u(b) = 0
    u_a = 0
    u_b = 0


    # The indexed variables
    # x_j = a + j * deltaX
    # t_m = m * deltaT


    # Initial Data u_0(x) = -10*(x-a)(x-b)
    u_initial = lambda x: -10 *(x-a)*(x-b)
    u = np.array([[u_a] + [u_initial(a + j * deltaX) for j in range(1,J)] + [u_b]], dtype=float)
    for m in range(M):
        ump1 = np.array([[u_a] + [0 for _ in range(1, J)] + [u_b]], dtype=float)
        # Explicit Scheme
        for j in range(1,J):
            ump1[0,j] = u[m,j] + cfl * (u[m,j+1] - 2 * u[m,j] + u[m,j-1])
        u = np.vstack([u, ump1])
    
    for m in range(0, M, 100):
        plt.plot([a + j * deltaX for j in range(J+1)], u[m])
    plt.show()

    print(u)


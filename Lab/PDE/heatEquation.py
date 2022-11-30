import numpy as np

if __name__ == "__main__":
    # Based on NSPDEs 5.5 Finite Difference Approximation of the Heat Equation
    # Parameters
    # \{(x,t) \in (a,b) \times [0,T]\}
    a = -1
    b = 1
    J = 20              # Final spacial point index
    T = 10              # Final time (in second)
    M = 10000           # Final time step
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    cfl = deltaT / (deltaX ** 2)            # Want this below 1/2
    # Dirichlet Boundary Conditions u(a) = 0, u(b) = 0
    u_a = 0
    u_b = 0


    # The indexed variables
    # x_j = a + j * deltaX
    # t_m = m * deltaT


    # Initial Data u_0(x) = -10*(x-a)(x-b)
    u_initial = lambda x: -10 *(x-a)*(x-b)
    u = np.array([[u_a] + [u_initial(a + j * deltaX) for j in range(1,J)] + [u_b]])
    for m in range(M):
        ump1 = np.array([[u_a] + [0 for _ in range(1, J)] + [u_b]])
        # Explicit Scheme
        for j in range(1,J):
            ump1[0,j] = u[m,j] + cfl * (u[m,j+1] - 2 * u[m,j] + u[m,j-1])
        u = np.vstack([u, ump1])

    print(u)


import matplotlib.pyplot as plt
from discrete import *
#import numpy as np

# Based on NSPDEs 5.5 Finite Difference Approximation of the Heat Equation
# Parameters
# \{(x,t) \in (a,b) \times [0,T]\}
# J Final spacial point index
# T Final time (in second)
# M Final time step
# deltaT Time mesh size
# deltaX Space mesh size
# Dirichlet Boundary Conditions u_a = u(a), u_b = u(b)

# This returns an array of u values at each time step and spacial position

# This returns an array of u values at each time step and spacial position (Periodic Condition)
def biharmonicEvolveExplicitEulerPeriodic(a=-1,b=1,J=20,T=10000,M=10000000,u_initial = lambda x: 10):
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    pseudocfl = deltaT / (deltaX ** 4)
    if(pseudocfl > 0.5):
        print(f"Warning: the CFL number is {pseudocfl}, so stability is not guaranteed")
    
    # Explicit Scheme
    u = [[u_initial(a + j * deltaX) for j in range(0,J + 1)]]
    for m in range(M):
        ump1 = [[0 for _ in range(0, J + 1)]]
        # Explicit Scheme
        for j in range(2,J-1):
            ump1[0][j] = u[m][j] - pseudocfl * (u[m][j+2] - 4 * u[m][j+1] + 6 * u[m][j] - 4 * u[m][j-1] + u[m][j-2])
        ump1[0][0] = u[m][0] - pseudocfl * (u[m][2] - 4 * u[m][1] + 6 * u[m][0] - 4 * u[m][J] + u[m][J-1])
        ump1[0][1] = u[m][1] - pseudocfl * (u[m][3] - 4 * u[m][2] + 6 * u[m][1] - 4 * u[m][0] + u[m][J])
        ump1[0][J-1] = u[m][J-1] - pseudocfl * (u[m][0] - 4 * u[m][J] + 6 * u[m][J-1] - 4 * u[m][J-2] + u[m][J-3])
        ump1[0][J] = u[m][J] - pseudocfl * (u[m][1] - 4 * u[m][0] + 6 * u[m][J] - 4 * u[m][J-1] + u[m][J-2])
        u += ump1
    return u

if __name__ == "__main__":
    a = -1
    b = 1
    J = 20                 # Final spacial point index
    T = 100                # Final time (in second)
    M = 10000000             # Final time step
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    u_initial = lambda x : -(x-1)*(x+1)*10
    u = biharmonicEvolveExplicitEulerPeriodic(a, b, J, T, M, u_initial=u_initial)
    deltaX = (b-a) / J  # Space mesh size
    
    plt.rcParams['text.usetex'] = True

    for m in range(0, M, 10):
        plt.plot([a + j * deltaX for j in range(J+1)], u[m])
    #axs[1].plot(range(M), [discretelPNorm(u[m], deltaX, p=1) for m in range(M)])
    #fig.title(r"$u_{tt} = \Delta u$ in $(x,t) \in \left[" + f"({a},{b}), (0,{T})" + r"\right]$")
    #fig.suptitle(f"{J + 1} spacial mesh points, {M+1} time mesh points, {J * (M+1)} points overall")
    plt.suptitle(r"$u_{tt} = -\Delta^{2} u$ in $(x,t) \in \left[" + f"({a},{b}), (0,{T})" + r"\right]$" + "\n" + f"{J + 1} spacial mesh points, {M+1} time mesh points, {J * (M+1)} points overall")
    plt.show()

    print(u[-1])


import matplotlib.pyplot as plt

# Parameters
# \{(x,t) \in (a,b) \times [0,T]\}
# J Final spacial point index
# T Final time (in second)
# M Final time step
# deltaT Time mesh size
# deltaX Space mesh size
def nonlinearEvolve1(a=-1,b=1,J=20,T=10000,M=10000000,u_initial = lambda x: 1 if (x < 1 and -1 < x) else 0):
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    cfl = deltaT / (deltaX ** 2)            # Want this below 1/2
    if(cfl > 0.5):
        print(f"Warning: the CFL number is {cfl}, so stability is not guaranteed")
    
    # Explicit Scheme
    u = [[u_initial(a + j * deltaX) for j in range(0,J + 1)]]
    for m in range(M):
        ump1 = [[0 for _ in range(0, J + 1)]]
        # Explicit Scheme
        for j in range(1,J):
            ump1[0][j] = u[m][j] + cfl * (2 * u[m][j] * (u[m][j+1] - 2 * u[m][j] + u[m][j-1]) + 0.5 * ((u[m][j+1] - u[m][j-1])**2))
        ump1[0][0] = u[m][j] + cfl * (2 * u[m][0] * (u[m][1] - 2 * u[m][0] + u[m][J]) + 0.5 * ((u[m][1] - u[m][J])**2))
        ump1[0][J] = u[m][j] + cfl * (2 * u[m][J] * (u[m][0] - 2 * u[m][J] + u[m][J-1]) + 0.5 * ((u[m][0] - u[m][J-1])**2))
        u += ump1
    return u

def nonlinearEvolve2(a=-1,b=1,J=20,T=10000,M=10000000,u_initial = lambda x: 1 if (x < 1 and -1 < x) else 0):
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    cfl = deltaT / (deltaX ** 2)            # Want this below 1/2
    if(cfl > 0.5):
        print(f"Warning: the CFL number is {cfl}, so stability is not guaranteed")
    
    # Explicit Scheme
    u = [[u_initial(a + j * deltaX) for j in range(0,J + 1)]]
    for m in range(M):
        ump1 = [[0 for _ in range(0, J + 1)]]
        # Explicit Scheme
        for j in range(1,J):
            #ump1[0][j] = u[m][j] + cfl * (2 * u[m][j] * (u[m][j+1] - 2 * u[m][j] + u[m][j-1]) + 0.5 * ((u[m][j+1] - u[m][j-1])**2))
            ump1[0][j] = u[m][j] + cfl * (u[m][j+1] ** 2 - 2 * u[m][j] ** 2 + u[m][j-1] ** 2)
        ump1[0][0] = u[m][0] + cfl * (u[m][1] ** 2 - 2 * u[m][0] ** 2 + u[m][J] ** 2)
        ump1[0][J] = u[m][J] + cfl * (u[m][0] ** 2 - 2 * u[m][J] ** 2 + u[m][J-1] ** 2)
        #ump1[0][0] = u[m][j] + cfl * (2 * u[m][0] * (u[m][1] - 2 * u[m][0] + u[m][J]) + 0.5 * ((u[m][1] - u[m][J])**2))
        #ump1[0][J] = u[m][j] + cfl * (2 * u[m][J] * (u[m][0] - 2 * u[m][J] + u[m][J-1]) + 0.5 * ((u[m][0] - u[m][J-1])**2))
        u += ump1
    return u

if __name__ == "__main__":
    a = -1
    b = 1
    J = 40                 # Final spacial point index
    T = 10                # Final time (in second)
    M = 80000             # Final time step
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size
    u_initial = lambda x: 1 if (x < 0.2 and -0.2 < x) else 0
    u = nonlinearEvolve2(a, b, J, T, M, u_initial)
    deltaX = (b-a) / J  # Space mesh size
    
    plt.rcParams['text.usetex'] = True

    for m in range(0, M, 10):
        plt.plot([a + j * deltaX for j in range(J+1)], u[m])
    plt.suptitle(r"$u_{t} = \Delta (u^2)$ in $(x,t) \in \left[" + f"({a},{b}), (0,{T})" + r"\right]$" + "\n" + f"{J + 1} spacial mesh points, {M+1} time mesh points, {J * (M+1)} points overall")
    plt.show()


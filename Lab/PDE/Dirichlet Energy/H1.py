import matplotlib.pyplot as plt

# u_t = -u
# J Final spacial point index
# T Final time (in second)
# M Final time step
# deltaT Time mesh size
# deltaX Space mesh size
# Dirichlet Boundary Conditions u_a = u(a), u_b = u(b)
def ode1EvolveExplicitEulerDirichlet(a=-1, b=1, J=20, T = 10000, M = 10000000, u_a=0, u_b=0, u_initial = lambda x: 10):
    deltaT = T/M        # Time mesh size
    deltaX = (b-a) / J  # Space mesh size

    # Explicit Scheme
    u = [[u_a] + [u_initial(a + j * deltaX) for j in range(1,J)] + [u_b]]
    for m in range(M):
        ump1 = [[u_a] + [0 for _ in range(1, J)] + [u_b]]
        # Explicit Scheme
        for j in range(1,J):
            ump1[0][j] = u[m][j] * (1 - deltaT)
        u += ump1
    return u

if __name__ == "__main__":
    a = -1
    b = 1
    J = 20              # Final spacial point index
    T = 10000           # Final time (in second)
    M = 10000000        # Final time step
    u_initial = lambda x : -10 * (x-a) * (x-b)
    u = ode1EvolveExplicitEulerDirichlet(a,b,J,T,M,u_a=0, u_b=0, u_initial=u_initial)
    deltaX = (b-a) / J  # Space mesh size


    plt.rcParams['text.usetex'] = True

    for m in range(0, M, 50):
        plt.plot([a + j * deltaX for j in range(J+1)], u[m])
    plt.title(r"$u_{tt} = u$ in $(x,t) \in \left[" + f"({a},{b}), (0,{T})" + r"\right]$")
    plt.suptitle(f"{J + 1} spacial mesh points, {M+1} time mesh points, {J * (M+1)} points overall")
    plt.show()


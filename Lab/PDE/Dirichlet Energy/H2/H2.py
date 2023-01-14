
# u_t u = \Delta^{-1} u
#
# Equivalently, writing \Delta (u_t) = u

# J Final spacial point index
# T Final time (in second)
# M Final time step
# deltaT Time mesh size
# deltaX Space mesh size
# Dirichlet Boundary Conditions u_a = u(a), u_b = u(b)
def inverseLaplacianFlowEuler(a=-1, b=1, J=20, T = 100, M = 100000, u_a=0, u_b=0, u_initial = lambda x: 10):
    pass

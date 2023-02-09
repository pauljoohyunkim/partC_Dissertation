import numpy as np
import matplotlib.pyplot as plt

# Steepest Descent Method
# f: Function (Taking a numpy array as its argument)
# df: Gradient of f (Again taking a numpy array as its argument), outputting a numpy array)
# alpha: Fixed step size
# x0: Initial argument (List)
# visualize_x: Python function for outputting the values. If not None, visualize_x(x) will be called at every iterate
# (Eg: visualize_x=print)
# visualize_f: Python function for outputting the function evaluated value. If not none,
# visualize_f(x) will be called every iteration.
# (Eg: visualize_f=print)
def sdm(f, df, alpha, x0, n=1000, visualize_x=None, visualize_f=None):
    x = x0
    for t in range(n):
        gradient = df(x)
        x = x - alpha * gradient
        if visualize_x:
            visualize_x(x)
        if visualize_f:
            visualize_f(f(x))



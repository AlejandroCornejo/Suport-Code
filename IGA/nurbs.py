
import numpy as np
import matplotlib.pyplot as pl

'''
In this file we create a set of Spline basis functions used in IGA
Then we use these functions to develop a rod IGA element

knot_vector = [x0, x1, ..., xi_{n+p+1}]

p: order of the polynomials
n: number of control points === number of basis functions
m: number of knots

m = n + p + 1

Ref: ISOGEOMETRIC STRUCTURAL ANALYSIS AND DESIGN, R. Wüchner, TUM
'''


def N0_BasisFunctions(xi, KnotVector):
    # returns the values of the shape functions for a give xi
    knot_size = KnotVector.size
    N0 = np.zeros(knot_size-1)

    for i in range(knot_size-1):
        if xi >= KnotVector[i] and xi < KnotVector[i + 1]:
            N0[i] = 1.0
        else:
            N0[i] = 0.0
    return N0

# local coordinates of the line
xi_vector = np.linspace(0.0, 1.0, 500)

# Let's try to reproduce Wüchner pg 18 Fig. 


# C0 basis functions
knot = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

# The number of basis functions is n = m - p - 1 = m - 1
N0 = np.zeros([knot.size-1, xi_vector.size]) # row if each N0_i, columns are values of xi


for counter, xi in enumerate(xi_vector):
    N0[:, counter] = N0_BasisFunctions(xi, knot)


# pl.plot(xi_vector, N0[0, :] , label="N0_0")
# pl.plot(xi_vector, N0[1, :] , label="N0_1")
# pl.plot(xi_vector, N0[2, :] , label="N0_2")
# pl.plot(xi_vector, N0[3, :] , label="N0_3")
pl.plot(xi_vector, N0[4, :] , label="N0_4")





pl.grid()
pl.legend()
pl.show()
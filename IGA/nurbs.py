
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

#----------------------------------------------------------------
def N0_BasisFunctions(xi, KnotVector):
    # returns the values of the shape functions for a give xi
    knot_size = KnotVector.size
    p = 0
    n = knot_size - 1 - p
    N0 = np.zeros(n)

    for i in range(n):
        if xi >= KnotVector[i] and xi < KnotVector[i + 1]:
            N0[i] = 1.0
        else:
            N0[i] = 0.0
    return N0
#----------------------------------------------------------------
def N1_BasisFunctions(xi, KnotVector):
    knot_size = KnotVector.size
    N0 = N0_BasisFunctions(xi, KnotVector)
    p = 1
    n = knot_size - 1 - p
    N1 = np.zeros(n)
    for i in range(n):
        N1[i] = (xi - KnotVector[i]) / (KnotVector[i+p]-KnotVector[i]) * N0[i] + \
                (KnotVector[i+p+1] - xi) / (KnotVector[i+p+1]-KnotVector[i+1]) * N0[i+1]
    return N1
#----------------------------------------------------------------
def N2_BasisFunctions(xi, KnotVector):
    knot_size = KnotVector.size
    N1 = N1_BasisFunctions(xi, KnotVector)
    p = 2
    n = knot_size - 1 - p
    N2 = np.zeros(n)
    for i in range(n):
        N2[i] = (xi - KnotVector[i]) / (KnotVector[i+p]-KnotVector[i]) * N1[i] + \
                (KnotVector[i+p+1] - xi) / (KnotVector[i+p+1]-KnotVector[i+1]) * N1[i+1]
    return N2
#----------------------------------------------------------------

# local coordinates of the line
xi_vector = np.linspace(0.0, 1.0, 2000)

# Let's try to reproduce Wüchner pg 18 Fig. 


# C0 basis functions
# knot = np.linspace(0.0, 1.0, 10)
knot = np.array([0.0, 0.1 , 0.12, 0.3, 0.4, 0.5, 0.52, 0.7, 0.8, 0.9, 1.0])

# The number of basis functions is n = m - p - 1 = m - 1
N0 = np.zeros([knot.size-1,   xi_vector.size]) # row if each N0_i, columns are values of xi
N1 = np.zeros([knot.size-1-1, xi_vector.size]) # row if each N1_i, columns are values of xi
N2 = np.zeros([knot.size-1-2, xi_vector.size]) # row if each N1_i, columns are values of xi


for counter, xi in enumerate(xi_vector):
    N0[:, counter] = N0_BasisFunctions(xi, knot)

for counter, xi in enumerate(xi_vector):
    N1[:, counter] = N1_BasisFunctions(xi, knot)

for counter, xi in enumerate(xi_vector):
    N2[:, counter] = N2_BasisFunctions(xi, knot)


# for i in range(knot.size - 1 - 1):
#     pl.plot(xi_vector, N1[i, :] , label="N_1_" + str(i))

for i in range(knot.size - 1 - 2):
    pl.plot(xi_vector, N2[i, :] , label="N_2_" + str(i))



pl.grid()
pl.legend()
pl.show()
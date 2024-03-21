
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
        denom1 = (KnotVector[i+p]-KnotVector[i])
        denom2 = (KnotVector[i+p+1]-KnotVector[i+1])
        if denom1 != 0.0:
            N1[i] += (xi - KnotVector[i]) / denom1 * N0[i]
        if denom2 != 0.0:
            N1[i] += (KnotVector[i+p+1] - xi) / denom2 * N0[i+1]
    return N1
#----------------------------------------------------------------
def N2_BasisFunctions(xi, KnotVector):
    knot_size = KnotVector.size
    N1 = N1_BasisFunctions(xi, KnotVector)
    p = 2
    n = knot_size - 1 - p
    N2 = np.zeros(n)
    for i in range(n):
        denom1 = (KnotVector[i+p]-KnotVector[i])
        denom2 = (KnotVector[i+p+1]-KnotVector[i+1])
        if denom1 != 0.0:
            N2[i] += (xi - KnotVector[i]) / denom1 * N1[i]
        if denom2 != 0.0:
            N2[i] += (KnotVector[i+p+1] - xi) / denom2 * N1[i+1]
    return N2
#----------------------------------------------------------------
def N2_DerivativesBasisFunctions(xi, KnotVector):
    knot_size = KnotVector.size
    N1 = N1_BasisFunctions(xi, KnotVector)
    p = 2
    n = knot_size - 1 - p
    dN2_dx = np.zeros(n)
    for i in range(n):
        denom1 = (KnotVector[i+p]-KnotVector[i])
        denom2 = (KnotVector[i+p+1]-KnotVector[i+1])
        if denom1 != 0.0:
            dN2_dx[i] += p / denom1 * N1[i]
        if denom2 != 0.0:
            dN2_dx[i] -= p / denom2 * N1[i+1]
    return dN2_dx
#----------------------------------------------------------------
def N3_BasisFunctions(xi, KnotVector):
    knot_size = KnotVector.size
    N2 = N2_BasisFunctions(xi, KnotVector)
    p = 3
    n = knot_size - 1 - p
    N3 = np.zeros(n)
    for i in range(n):
        denom1 = (KnotVector[i+p]-KnotVector[i])
        denom2 = (KnotVector[i+p+1]-KnotVector[i+1])
        if denom1 != 0.0:
            N3[i] += (xi - KnotVector[i]) / denom1 * N2[i]
        if denom2 != 0.0:
            N3[i] += (KnotVector[i+p+1] - xi) / denom2 * N2[i+1]
    return N3
#----------------------------------------------------------------
def N4_BasisFunctions(xi, KnotVector):
    knot_size = KnotVector.size
    N3 = N3_BasisFunctions(xi, KnotVector)
    p = 4
    n = knot_size - 1 - p
    N4 = np.zeros(n)
    for i in range(n):
        denom1 = (KnotVector[i+p]-KnotVector[i])
        denom2 = (KnotVector[i+p+1]-KnotVector[i+1])
        if denom1 != 0.0:
            N4[i] += (xi - KnotVector[i]) / denom1 * N3[i]
        if denom2 != 0.0:
            N4[i] += (KnotVector[i+p+1] - xi) / denom2 * N3[i+1]
    return N4
#----------------------------------------------------------------

# local coordinates of the line
xi_vector = np.linspace(-1.0, 1.0, 2000)


# C0 basis functions
# knot = np.linspace(0.0, 1.0, 10)
# knot = np.array([0.0, 0.0 , 0.0, 0.0, 0.4, 0.5, 0.6, 1.0, 1.0, 1.0, 1.0])
knot = np.array([0, 0, 0,  2.0, 2.0, 2.0])

# The number of basis functions is n = m - p - 1 = m - 1
N0 = np.zeros([knot.size-1,   xi_vector.size]) # row if each N0_i, columns are values of xi
N1 = np.zeros([knot.size-1-1, xi_vector.size]) # row if each N1_i, columns are values of xi
N2 = np.zeros([knot.size-1-2, xi_vector.size]) # row if each N1_i, columns are values of xi
dN2_dx = np.zeros([knot.size-1-2, xi_vector.size]) # row if each N1_i, columns are values of xi
N3 = np.zeros([knot.size-1-3, xi_vector.size]) # row if each N1_i, columns are values of xi
N4 = np.zeros([knot.size-1-4, xi_vector.size]) # row if each N1_i, columns are values of xi


# for counter, xi in enumerate(xi_vector):
#     N0[:, counter] = N0_BasisFunctions(xi, knot)
# for i in range(knot.size - 1 - 0):
#     pl.plot(xi_vector, N0[i, :] , label="N_0_" + str(i))


# for counter, xi in enumerate(xi_vector):
#     N1[:, counter] = N1_BasisFunctions(xi, knot)
# for i in range(knot.size - 1 - 1):
#     pl.plot(xi_vector, N1[i, :] , label="N_1_" + str(i))


# for counter, xi in enumerate(xi_vector):
#     N2[:, counter] = N2_BasisFunctions(xi, knot)
# for counter, xi in enumerate(xi_vector):
#     dN2_dx[:, counter] = N2_DerivativesBasisFunctions(xi, knot)

# for i in range(knot.size - 1 - 2):
#     pl.plot(xi_vector, N2[i, :] , label="N_2_" + str(i))

# for i in range(knot.size - 1 - 2):
#     pl.plot(xi_vector, dN2_dx[i, :] , label="dN_2_" + str(i))


# for counter, xi in enumerate(xi_vector):
#     N3[:, counter] = N3_BasisFunctions(xi, knot)
# for i in range(knot.size - 1 - 3):
#     pl.plot(xi_vector, N3[i, :] , label="N_3_" + str(i))

# for counter, xi in enumerate(xi_vector):
#     N4[:, counter] = N4_BasisFunctions(xi, knot)
# for i in range(knot.size - 1 - 4):
#     pl.plot(xi_vector, N4[i, :] , label="N_4_" + str(i))

# pl.grid()
# pl.legend()
# pl.show()


'''
Here we begin to compute the K of a 3 Control Points IGA patch with
quadratic basis functions, 2 elements,  knots = [0,0,0,1,1,1]
      0           1           2
    O-----------O-----------O

The patch is separated into 2 "elements" inside of which the
u = Sum(Ni * ui) for i = 0,...,n

Each element is integrated using Gauss quadratures locally

'''

# The local axes xi is [-1, 1]
xi_ip   = 0.0
weight = 2.0

E = 2.1e9
A = 0.001
L = 2.0


# Element 1
L1 = 1.0
xc1 = 0.5
J1 = 0.5*L1
x_ip = L1*xi_ip*0.5 + xc1 # relation between xi and x, xi = 2 / L * (x - xc)
K1 = np.zeros([3, 3])
N2_derivatives = N2_DerivativesBasisFunctions(x_ip, knot) # These derivatives are in the parametric space
K1 += E * A * np.outer(N2_derivatives, N2_derivatives) * J1 * weight # Jacobian from the parametric to the natural coords

# Element 2
L2 = 1.0
xc2 = 1.5
J2 = 0.5*L2
x_ip = L2*xi_ip*0.5 + xc2 # relation between xi and x
K2 = np.zeros([3, 3])
N2_derivatives = N2_DerivativesBasisFunctions(x_ip, knot)
K2 += E * A * np.outer(N2_derivatives, N2_derivatives) * J2 * weight

# The assembly is trivial
K = K1+K2

K_bc = np.zeros([2,2])

K_bc[0,0] = K[1,1]
K_bc[1,1] = K[2,2]
K_bc[0,1] = K[1, 2]
K_bc[1,0] = K[2, 1]

p = np.array([0.0, 1e3])

displ = np.dot(np.linalg.inv(K_bc), p)

print(displ)

analytical_displ = p[1] * L / A / E
print(analytical_displ)
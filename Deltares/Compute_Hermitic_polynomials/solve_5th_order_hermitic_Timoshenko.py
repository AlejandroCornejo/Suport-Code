import sympy as sp
import numpy as np

'''
In this script I am solving the required system of equations to obtain the alpha parameters
such that:
    v = alpha_0 + alpha_1*xi + alpha_2*xi**2 + alpha_3*xi**3 + alpha_4*xi**4 + alpha_5*xi**5

fulfiling the boundary conditions at the 3 nodes

Timoshenko theory: theta = dv/dx + gamma

phi is the shear parameter
'''

v1, theta1, v2, theta2, v3, theta3, L, a0, a1, a2, a3, a4, a5, phi, xi = sp.symbols('v1 theta1 v2 theta2 v3 theta3 L a0 a1 a2 a3 a4 a5 phi xi')

eq1 = sp.Eq(a0 - a1 + a2 - a3 + a4 - a5, v1)
eq2 = sp.Eq(a0                         , v2)
eq3 = sp.Eq(a0 + a1 + a2 + a3 + a4 + a5, v3)

eq4 = sp.Eq(2.0 / L * (a1 -2.0 * a2 + 3.0 * a3 - 4.0 * a4 + 5.0 * a5) + phi * L**2 / 12.0 * 8.0 / L**3 * (6.0 * a3 - 24.0 * a4 + 60.0 * a5), theta1) # v' + gamma
eq5 = sp.Eq(2.0 / L * (a1)                                            + phi * L**2 / 12.0 * 8.0 / L**3 * (6.0 * a3)                        , theta2) # v' + gamma
eq6 = sp.Eq(2.0 / L * (a1 +2.0 * a2 + 3.0 * a3 + 4.0 * a4 + 5.0 * a5) + phi * L**2 / 12.0 * 8.0 / L**3 * (6.0 * a3 + 24.0 * a4 + 60.0 * a5), theta3) # v' + gamma

solution = sp.solve((eq1, eq2, eq3, eq4, eq5, eq6), (a0, a1, a2, a3, a4, a5))

# print("The value of alpha_0 is: ", solution[a0])
# print("The value of alpha_1 is: ", solution[a1])
# print("The value of alpha_2 is: ", solution[a2])
# print("The value of alpha_3 is: ", solution[a3])
# print("The value of alpha_4 is: ", solution[a4])
# print("The value of alpha_5 is: ", solution[a5])

'''
Output:
    The value of alpha_0 is:  v2
    The value of alpha_1 is:  (-L*phi*theta1 - 18.0*L*phi*theta2 - L*phi*theta3 - 2.0*L*theta2 - 40.0*phi**2*v1 + 40.0*phi**2*v3 - 10.0*phi*v1 + 10.0*phi*v3)/(80.0*phi**2 - 20.0*phi - 4.0)
    The value of alpha_2 is:  (L*theta1 - L*theta3 + 16.0*phi*v1 - 32.0*phi*v2 + 16.0*phi*v3 + 8.0*v1 - 16.0*v2 + 8.0*v3)/(32.0*phi + 8.0)
    The value of alpha_3 is:  (40.0*L*phi*theta2 + L*theta1 + 8.0*L*theta2 + L*theta3 + 40.0*phi*v1 - 40.0*phi*v3 + 10.0*v1 - 10.0*v3)/(160.0*phi**2 - 40.0*phi - 8.0)
    The value of alpha_4 is:  (-L*theta1 + L*theta3 - 4.0*v1 + 8.0*v2 - 4.0*v3)/(32.0*phi + 8.0)
    The value of alpha_5 is:  (2.0*L*phi*theta1 - 4.0*L*phi*theta2 + 2.0*L*phi*theta3 - L*theta1 - 4.0*L*theta2 - L*theta3 - 6.0*v1 + 6.0*v3)/(160.0*phi**2 - 40.0*phi - 8.0)
'''


# let's check...
L = 1.0
phi = 0.0
# v1     = 0
# v2     = 1
# v3     = 0.0
# theta1 = 0
# theta2 = 0
# theta3 = 0

# a_0 = v2
# a_1 = (-L*phi*theta1 - 18.0*L*phi*theta2 - L*phi*theta3 - 2.0*L*theta2 - 40.0*phi**2*v1 + 40.0*phi**2*v3 - 10.0*phi*v1 + 10.0*phi*v3)/(80.0*phi**2 - 20.0*phi - 4.0)
# a_2 =  (L*theta1 - L*theta3 + 16.0*phi*v1 - 32.0*phi*v2 + 16.0*phi*v3 + 8.0*v1 - 16.0*v2 + 8.0*v3)/(32.0*phi + 8.0)
# a_3 = (40.0*L*phi*theta2 + L*theta1 + 8.0*L*theta2 + L*theta3 + 40.0*phi*v1 - 40.0*phi*v3 + 10.0*v1 - 10.0*v3)/(160.0*phi**2 - 40.0*phi - 8.0)
# a_4 = (-L*theta1 + L*theta3 - 4.0*v1 + 8.0*v2 - 4.0*v3)/(32.0*phi + 8.0)
# a_5 = (2.0*L*phi*theta1 - 4.0*L*phi*theta2 + 2.0*L*phi*theta3 - L*theta1 - 4.0*L*theta2 - L*theta3 - 6.0*v1 + 6.0*v3)/(160.0*phi**2 - 40.0*phi - 8.0)

'''
Let's collect the indivual polynomials N = [N1, N1_bar, N2, N2_bar, N3, N3_bar]
    v == [a0, a1, a2, a3, a4, a5] * [1, xi, xi**2, xi**3, xi**4, xi**5]

    We want to transform it to v == [N1, N1_bar, N2, N2_bar, N3, N3_bar] * [v1, theta1, v2, theta2, v3, theta3]

-> v1
N1 = (-40.0*phi**2 - 10.0*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + 
     (16.0*phi + 8.0 )/(32.0*phi + 8.0) * xi**2 + 
     (40.0*phi + 10.0 )/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 +
     (- 4.0 )/(32.0*phi + 8.0) * xi**4 +
     (- 6.0 )/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

-> theta1
N1_bar =    (-L*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + 
            (L)/(32.0*phi + 8.0) * xi**2 + 
            (L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + 
            (-L)/(32.0*phi + 8.0) * xi**4 + 
            (2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

-> v2
N2 =  1.0 + 
      (- 32.0*phi  - 16.0 )/(32.0*phi + 8.0) * xi**2 +
      (8.0)/(32.0*phi + 8.0) * xi**4

-> theta2
N2_bar = (- 18.0*L*phi - 2.0*L ) / (80.0*phi**2 - 20.0*phi - 4.0) * xi +
         (40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 +
         (- 4.0*L*phi  - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

-> v3
N3 = (40.0*phi**2 + 10.0*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi +
     (16.0*phi + 8.0)/(32.0*phi + 8.0) * xi**2 +
     (- 40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 +
     (- 4.0)/(32.0*phi + 8.0) * xi**4 +
     (6.0)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

-> theta3
N3_bar = (- L*phi ) / (80.0*phi**2 - 20.0*phi - 4.0) * xi +
         (- L)/(32.0*phi + 8.0) * xi**2 +
         (L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + 
         (L)/(32.0*phi + 8.0) * xi**4 +
         (2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5
'''


# def ComputeV(a0, a1, a2, a3, a4, a5, L, xi):
#     return a0 + a1*xi + a2*xi**2 + a3*xi**3 + a4*xi**4 + a5*xi**5

# def ComputeTheta(a0, a1, a2, a3, a4, a5, L, phi, xi):
#     return 2.0 / L * (a1 + 2.0 * a2*xi + 3.0*a3*xi**2 + 4.0*a4*xi**3 + 5.0*a5*xi**4) + phi * L**2 / 12.0 * 8 / L**3 * (6*a3 + 24*a4*xi+60*a5*xi**2)

# print("\n")
# print("The calculated v at xi=-1 is: ", ComputeV(a_0, a_1, a_2, a_3, a_4, a_5, L, -1))
# print("The calculated v at xi= 1 is: ", ComputeV(a_0, a_1, a_2, a_3, a_4, a_5, L,  1))
# print("The calculated v at xi= 0 is: ", ComputeV(a_0, a_1, a_2, a_3, a_4, a_5, L,  0))
# print("The calculated v at xi= 1 is: ", ComputeV(alpha_0, alpha_1, alpha_2, alpha_3,  1))
# print("The calculated Theta at xi=-1 is: ", ComputeTheta(a_0, a_1, a_2, a_3, a_4, a_5, L, phi, -1))

def ComputeN(L, phi, xi):
    N = np.zeros(6)
    N[0] = (-40.0*phi**2 - 10.0*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
            (16.0*phi + 8.0 )/(32.0*phi + 8.0) * xi**2 + \
            (40.0*phi + 10.0 )/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
            (- 4.0 )/(32.0*phi + 8.0) * xi**4 + \
            (- 6.0 )/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

    N[1] = (-L*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
            (L)/(32.0*phi + 8.0) * xi**2 + \
            (L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
            (-L)/(32.0*phi + 8.0) * xi**4 + \
            (2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5
    N[2] =  1.0 + \
      (- 32.0*phi  - 16.0 )/(32.0*phi + 8.0) * xi**2 + \
      (8.0)/(32.0*phi + 8.0) * xi**4

    N[3] = (- 18.0*L*phi - 2.0*L ) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
         (40.0*L*phi + 8.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
         (- 4.0*L*phi  - 4.0*L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

    N[4] = (40.0*phi**2 + 10.0*phi) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
            (16.0*phi + 8.0)/(32.0*phi + 8.0) * xi**2 + \
            (- 40.0*phi - 10.0)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 +\
            (- 4.0)/(32.0*phi + 8.0) * xi**4 +\
            (6.0)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5

    N[5] = (- L*phi ) / (80.0*phi**2 - 20.0*phi - 4.0) * xi + \
         (- L)/(32.0*phi + 8.0) * xi**2 + \
         (L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**3 + \
         (L)/(32.0*phi + 8.0) * xi**4 + \
         (2.0*L*phi - L)/(160.0*phi**2 - 40.0*phi - 8.0) * xi**5
    return N

# [v1, theta1, v2, theta2, v3, theta3]
u_e = np.zeros(6)
u_e[0] = 0.125
u_e[2] = -0.5
print("The calculated v at xi=-1 is: ", np.dot(ComputeN(L, phi, -1), u_e))
print("The calculated v at xi= 0 is: ", np.dot(ComputeN(L, phi,  0), u_e))

import numpy as np
import matplotlib.pyplot as pl
import math
from timoshenko_element import*

"""

Felippa and Oñate, "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability", Archives of Comp. Methods in Eng. (2021) 28:2021-2080

In this file we derive the main equations of 2-noded Timoshenko FE, firstly using a closed stiffness form (K) and then by integration (int BtDB dV)

"""

# We start the calculations
"""
    ^ v1       ^ v2
    |          |
    O----------O --> x , rot in 1 and 2

    u^(e) = {v1, rot1, v2, rot2}
"""

# Material properties
E = 2.1e9  # Pa
nu = 0.3   # -
A = 0.001  # m2
I = 1.0e-4 # m4

P = -10e3 # N

# Nodes
node_1 = Node(0.0, 0.0)
node_2 = Node(10.0, 0.0)

element = TimoshenkoElement2D2N(node_1, node_2, E, I, nu, A)

# K = element.CalculateDirectK()

# v1 = rot1 = 0
# K_bc = K[2:4, 2:4]

# f = np.array([P, 0.0]) # vertical descending force in node 2


# Solve the system
# u_bc = np.dot(np.linalg.inv(K_bc), f)

# print(u_bc)
# print("v2   = ", u_bc[0, 0])
# print("rot2 = ", u_bc[0, 1])

# analytical_v2 = P*element.Length * (1.0 / (element.G*element.As) + (3*element.Length**2-element.Length**2) / (6.0*E*I))

# print("Analytical v2 = ", analytical_v2)



u_e_test = np.zeros(4)
u_e_test[0] = 0.0  # v1
u_e_test[1] = 1.0  # rot1
u_e_test[2] = 0.0  # v2
u_e_test[3] = 1.0  # rot3

# element.PrintShapeFunctions()

# element.PrintDeflectionCurveFromNodalValues(u_e_test)
# element.PrintRotationCurveFromNodalValues(u_e_test)
# element.PrintStrainKinematics(u_e_test)
# element.Phi = 0.0

K_integrated = element.CalculateStiffnessMatrix(3)
# K = element.CalculateDirectK()

# print(K)
print(K_integrated)



import numpy as np
import matplotlib.pyplot as pl
import math
import builder_and_solver as BaS

from timoshenko_element_2D_2N import*
from node_2d import*

"""
Felippa and Oñate, "Accurate Timoshenko Beam Elements For Linear Elastostatics and LPB Stability", Archives of Comp. Methods in Eng. (2021) 28:2021-2080
In this file we derive the main equations of 2-noded Timoshenko FE, firstly using a closed stiffness form (K) and then by integration (int BtDB dV)
"""

# Material properties
E = 2.1e9  # Pa
nu = 0.3   # -
A = 116.0e-4  # m2
I = 48200.0e-8 # m4

# Loads
P = -1e3 # N
M = 5e3   # Nm

# Nodes
node_1 = Node2D(1, 0.0, 0.0)
node_2 = Node2D(2, 5.0, 0.0)

# Element
element = TimoshenkoElement2D2N(node_1, node_2, E, I, nu, A)
K_integrated = element.CalculateStiffnessMatrix(2)


############################################################
# NUMERICAL EXAMPLES
############################################################

# Case a) cantilever with axial P
# external force vector f = {Px1, Py1, M1, Px2, Py2, M2}

f = np.array([0.0, 0.0, 0.0, P, 0.0, 0.0])
fixed_dofs = [0,1,2]
K_bc = element.ApplyBoundaryConditionsToK(K_integrated, fixed_dofs)
displacement_vector = BaS.SolveSystem(K_bc, f)
analytical_displ = f[3] / A / E * element.Length
print("\n--> Case a), Cantilever with a horizontal/axial load")
print("The analytical axial displacement is: ", '{:.5e}'.format(analytical_displ), "m")
print("The FEM        axial displacement is: ", '{:.5e}'.format(displacement_vector[3]), "m") # OK

reactions = np.dot(K_integrated, displacement_vector)
print("The axial reaction is: ", '{:.5e}'.format(reactions[0]), " N\n")

# Case b) cantilever with vertical P
# external force vector f = {Px1, Py1, M1, Px2, Py2, M2}

f = np.array([0.0, 0.0, 0.0, 0.0, P, 0.0])
K_bc = element.ApplyBoundaryConditionsToK(K_integrated, [0,1,2])
displacement_vector = BaS.SolveSystem(K_bc, f)
analytical_displ = P*element.Length * (1.0 / (element.G*element.As) + (3*element.Length**2-element.Length**2) / (6.0*E*I))
print("--> Case b), Cantilever with a vertical load")
print("The analytical vert displacement is: ", '{:.5e}'.format(analytical_displ), "m")
print("The FEM        vert displacement is: ", '{:.5e}'.format(displacement_vector[4]), "m") # OK

reactions = np.dot(K_integrated, displacement_vector)
print("The axial reaction    is: ", '{:.5e}'.format(reactions[0]), " N")
print("The vertical reaction is: ", '{:.5e}'.format(reactions[1]), " N")
print("The moment reaction   is: ", '{:.5e}'.format(reactions[2]), " Nm\n")

# u_e = np.zeros(4)
# u_e[0] = 0.
# u_e[1] = 0.
# u_e[2] = displacement_vector[4]
# u_e[3] = displacement_vector[5]

# element.PrintDeflectionCurveFromNodalValues(u_e)
# element.PrintRotationCurveFromNodalValues(u_e)
# element.PrintBendingMomentsFromNodalValues(u_e)
# element.PrintStrainKinematics(u_e)
# element.PrintShearForceFromNodalValues(u_e)



# Case c) cantilever with axial and vertical load P, we fix the other extreme
# external force vector f = {Px1, Py1, M1, Px2, Py2, M2}

f = np.array([2*P, P, 0.0, 0.0, 0.0, 0.0])
K_bc = element.ApplyBoundaryConditionsToK(K_integrated, [3,4,5])
displacement_vector = BaS.SolveSystem(K_bc, f)
analytical_hor_displ = 2*P / A / E * element.Length
analytical_displ_vert = P*element.Length * (1.0 / (element.G*element.As) + (3*element.Length**2-element.Length**2) / (6.0*E*I))
print("--> Case c), Cantilever with a vertical and axial load")
print("The analytical hori displacement is: ", '{:.5e}'.format(analytical_hor_displ), "m")
print("The FEM        hori displacement is: ", '{:.5e}'.format(displacement_vector[0]), "m")
print("The analytical vert displacement is: ", '{:.5e}'.format(analytical_displ_vert), "m")
print("The FEM        vert displacement is: ", '{:.5e}'.format(displacement_vector[1]), "m")

reactions = np.dot(K_integrated, displacement_vector)
print("The axial reaction    is: ", '{:.5e}'.format(reactions[3]), " N")
print("The vertical reaction is: ", '{:.5e}'.format(reactions[4]), " N")
print("The moment reaction   is: ", '{:.5e}'.format(reactions[5]), " Nm\n")

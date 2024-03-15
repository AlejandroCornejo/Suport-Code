
import numpy as np
import matplotlib.pyplot as pl
import math
import builder_and_solver as BaS

from timoshenko_element_2D_2N import*
from timoshenko_element_2D_3N import*
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
L = 5.0

# Loads
P = -1e3 # N
M = 5e3   # Nm

# Nodes
node_1 = Node2D(1, 0.0,   0.0)
node_2 = Node2D(2, L,     0.0)
# node_2 = Node2D(2, 2.0*L,     0.0)
node_3 = Node2D(3, 2.0*L, 0.0)

# Element
element = TimoshenkoElement2D2N(node_1, node_2, E, I, nu, A)
K_integrated = element.CalculateStiffnessMatrix(2)




############################################################
# NUMERICAL EXAMPLES
############################################################

# Case a) cantilever with axial P
# external force vector f = {Px1, Py1, M1, Px2, Py2, M2}

# f = np.array([0.0, 0.0, 0.0, P, 0.0, 0.0])
# fixed_dofs = [0,1,2]
# K_bc = element.ApplyBoundaryConditionsToK(K_integrated, fixed_dofs)
# displacement_vector = BaS.SolveSystem(K_bc, f)
# analytical_displ = f[3] / A / E * element.Length
# print("\n--> Case a), Cantilever with a horizontal/axial load")
# print("The analytical axial displacement is: ", '{:.5e}'.format(analytical_displ), "m")
# print("The FEM        axial displacement is: ", '{:.5e}'.format(displacement_vector[3]), "m") # OK

# reactions = np.dot(K_integrated, displacement_vector)
# print("The axial reaction is: ", '{:.5e}'.format(reactions[0]), " N\n")




# Case b) cantilever with vertical P
# external force vector f = {Px1, Py1, M1, Px2, Py2, M2}

# f = np.array([0.0, 0.0, 0.0, 0.0, P, 0.0])
# K_bc = element.ApplyBoundaryConditionsToK(K_integrated, [0,1,2])
# displacement_vector = BaS.SolveSystem(K_bc, f)
# analytical_displ = P*element.Length * (1.0 / (element.G*element.As) + (3.*element.Length**2.-element.Length**2.) / (6.0*E*I))
# print("--> Case b), Cantilever with a vertical load")
# print("The analytical vert displacement is: ", '{:.5e}'.format(analytical_displ), "m")
# print("The FEM        vert displacement is: ", '{:.5e}'.format(displacement_vector[4]), "m") # OK

# reactions = BaS.ComputeReactions(K_integrated, displacement_vector)
# print("The axial reaction    is: ", '{:.5e}'.format(reactions[0]), " N")
# print("The vertical reaction is: ", '{:.5e}'.format(reactions[1]), " N")
# print("The moment reaction   is: ", '{:.5e}'.format(reactions[2]), " Nm\n")


# element.PrintDeflectionCurveFromNodalValues(displacement_vector)
# element.PrintRotationCurveFromNodalValues(displacement_vector)
# element.PrintBendingMomentsFromNodalValues(displacement_vector)
# element.PrintStrainKinematics(displacement_vector)
# element.PrintShearForceFromNodalValues(displacement_vector)

# Fint = element.CalculateInternalForcesVector(displacement_vector)
# Fext = np.dot(K_integrated, displacement_vector)
# print("The residual norm is: ", '{:.5e}'.format(np.linalg.norm(Fint-Fext)))



# Case c) cantilever with axial and vertical load P, we fix the other extreme
# external force vector f = {Px1, Py1, M1, Px2, Py2, M2}

# f = np.array([2*P, P, 0.0, 0.0, 0.0, 0.0])
# K_bc = element.ApplyBoundaryConditionsToK(K_integrated, [3,4,5])
# displacement_vector = BaS.SolveSystem(K_bc, f)
# analytical_hor_displ = 2*P / A / E * element.Length
# analytical_displ_vert = P*element.Length * (1.0 / (element.G*element.As) + (3*element.Length**2-element.Length**2) / (6.0*E*I))
# print("--> Case c), Cantilever with a vertical and axial load")
# print("The analytical hori displacement is: ", '{:.5e}'.format(analytical_hor_displ), "m")
# print("The FEM        hori displacement is: ", '{:.5e}'.format(displacement_vector[0]), "m")
# print("The analytical vert displacement is: ", '{:.5e}'.format(analytical_displ_vert), "m")
# print("The FEM        vert displacement is: ", '{:.5e}'.format(displacement_vector[1]), "m")

# reactions = np.dot(K_integrated, displacement_vector)
# print("The axial reaction    is: ", '{:.5e}'.format(reactions[3]), " N")
# print("The vertical reaction is: ", '{:.5e}'.format(reactions[4]), " N")
# print("The moment reaction   is: ", '{:.5e}'.format(reactions[5]), " Nm\n")





# element.RotateK(K_integrated)
# print(K_integrated)
# f = np.array([0.0, 0.0, 0.0, P, 0.0, 0.0])
# K_bc = element.ApplyBoundaryConditionsToK(K_integrated, [0,1,2])
# displacement_vector = BaS.SolveSystem(K_bc, f)
# print(displacement_vector)
# analytical_displ = P*element.Length * (1.0 / (element.G*element.As) + (3.*element.Length**2.-element.Length**2.) / (6.0*E*I))
# print("--> Case d), Rotated Cantilever with a vertical load")
# print("The analytical vert displacement is: ", '{:.5e}'.format(analytical_displ), "m")
# print("The FEM        vert displacement is: ", '{:.5e}'.format(displacement_vector[3]), "m") # OK






# 3 noded element

element_3N = TimoshenkoElement2D3N(node_1, node_2, node_3, E, I, nu, A)

K_3N = element_3N.CalculateStiffnessMatrix(4)
f = np.zeros(9)
f[7] = P
K_3N_bc = element.ApplyBoundaryConditionsToK(K_3N, [0,1,2])
displacement = BaS.SolveSystem(K_3N_bc, f)
# print(displacement)

# analytical_displ = P*element_3N.Length * ((3.*element_3N.Length**2.-element_3N.Length**2.) / (6.0*E*I)) # EB ok
analytical_displ = P*element_3N.Length * (1.0 / (element_3N.G*element_3N.As) + (3.*element_3N.Length**2.-element_3N.Length**2.) / (6.0*E*I))
print("--> Case b), Cantilever with a vertical load")
print("The analytical vert displacement is: ", '{:.5e}'.format(analytical_displ), "m")
print("The FEM        vert displacement is: ", '{:.5e}'.format(displacement[7]), "m")



# print(K_3N)
import Node
import Element
import Material
import Integration
import BuilderAndSolver as BaS
import BoundaryConditions as BC
import ReadFile
import scipy as sp 
import numpy as np
np.set_printoptions(precision = 3)

def Wait():
	input("Press Something")

FileName = input(" Name of the input file (with extension):   ")

if FileName == '0':
	FileName = "Input.dat"

File     = open(FileName, 'r')
ReadFile = ReadFile.ReadFile(File)

# Read materials, nodes, elements, FixedNodes and NodalForces to consider
Materials  = ReadFile.GetMaterials()
Nodes      = ReadFile.GetNodes()
Elements   = ReadFile.GetElements(Nodes, Materials)
FixedDoF   = ReadFile.GetFixedDoF(Nodes)
NodalForce = ReadFile.GetNodalForces(Nodes)

# General data
number_of_nodes    = ReadFile.GetNumberOfNodes()
number_of_elements = ReadFile.GetNumberOfElements()

# Builder and Solver
BuilderAndSolver = BaS.BuilderAndSolver(number_of_nodes)

# Integration Parameters
Integration = Integration.IntegrationQuadrature()
Integration.SetIntegrationPoints()
IntegrationPoints  = Integration.GetIntegrationPoints()
IntegrationWeigths = Integration.GetIntegrationWeights()

# Initialize Global LHS and RHS
GlobalLHS = np.zeros((number_of_nodes*3, number_of_nodes*3))
GlobalRHS = np.zeros(number_of_nodes*3)

# LOOP over elements to compute K in each element
for elem in range(number_of_elements):
	LocalLHS = np.zeros((24,24))
	# LOOP over Integration Points
	for IntPoint in range(Integration.GetNumberOfIP()):
		W  = IntegrationWeigths[IntPoint]
		E  = Elements[elem].GetMaterial().GetYoungModulus()
		nu = Elements[elem].GetMaterial().GetPoisson()
		C  = Elements[elem].GetLinearConstitutiveMatrix(E, nu)
		CoordIP = IntegrationPoints[IntPoint].GetCoordinates()

		N     = Elements[elem].GetLocalShapeFunctions(CoordIP)
		DN_De = Elements[elem].GetLocalShapeFunctionsGradients(CoordIP)
		J     = Elements[elem].GetJacobianMatrix(IntegrationPoints[IntPoint])
		detJ  = np.linalg.det(J)
		DN_DX = Elements[elem].CalculateDN_DX(DN_De, J)
		B     = Elements[elem].CalculateDeformationMatrix(DN_DX)
		Elements[elem].SetDeformationMatrix(B, IntPoint)

		LocalLHS += np.dot(np.dot(np.transpose(B),C),B)*W*detJ   # K = sum( Bt*C*B*W*det(J))
	BuilderAndSolver.AssembleLHS(GlobalLHS, LocalLHS, Elements[elem].GetDoF())

# The not-modified LHS is OriginalGLobalLHS
OriginalGlobalLHS = LocalLHS

# Apply Boundary Conditions
ImposedDisplacement = BC.BoundaryCondition("ImposedDisplacement", FixedDoF)
ImposedDisplacement.ApplyBoundaryCondition(GlobalLHS, GlobalRHS)

# Apply nodal forces
NodalForce = BC.BoundaryCondition("NodalForce", NodalForce)
NodalForce.ApplyBoundaryCondition(GlobalLHS, GlobalRHS)

# Solve the system
DisplVector = BuilderAndSolver.Solve(GlobalLHS, GlobalRHS)

# Assign the Displacements obtained to each Node
BuilderAndSolver.AssignDisplacementsToNodes(Nodes, DisplVector)

# Compute nodal Reactions and assign to nodes NO FUNCIONA BN
BuilderAndSolver.CalculateAndAssignReactions(OriginalGlobalLHS, DisplVector, Nodes)

# LOOP over elements to compute Strains and Stresses in each IP
for elem in range(number_of_elements):
	LocalLHS = np.zeros((24,24))
	# LOOP over Integration Points
	for IntPoint in range(Integration.GetNumberOfIP()):
		B = Elements[elem].GetDeformationMatrix(IntPoint)
		StrainVector = Elements[elem].CalculateStrainVector(B)
		Elements[elem].SetStrainVector(StrainVector, IntPoint)

		E  = Elements[elem].GetMaterial().GetYoungModulus()
		nu = Elements[elem].GetMaterial().GetPoisson()
		C  = Elements[elem].GetLinearConstitutiveMatrix(E, nu)

		StressVector = Elements[elem].CalculateStressVector(StrainVector, C)
		Elements[elem].SetStressVector(StressVector, IntPoint)

File.close()
# See results
#print(Nodes[0].GetReactions())
#print("Displ ",DisplVector)
print(Elements[0].GetStressVector(1))
#print(Elements[0].GetStressVector(1))
#print(Elements[elem].GetStrainVectors())

# Write Results TODO
OutFile = open ("Output.post.res", "w")
OutFile.write("RESULTADOS DEL ANALISIS \n")
for node in Nodes:
	OutFile.write(str(node.GetId())) 
	OutFile.write("   ")
	OutFile.write(str(node.GetX()))
	OutFile.write("   ")
	OutFile.write("   ")
	OutFile.write(str(node.GetY()))
	OutFile.write("   ")
	OutFile.write("   ")
	OutFile.write(str(node.GetZ()))
	OutFile.write("   ")
	OutFile.write("   \n")
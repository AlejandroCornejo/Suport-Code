import scipy as sp 
import numpy as np

def Wait():
	input("Press Something")

class BuilderAndSolver:
	def __init__(self, NumberOfNodes):
		self.__TotalDoF = 3 * NumberOfNodes
		self.__NumberOfNodes = NumberOfNodes

	def GetInfo(self):
		print("Builder and Solver")

	def AssembleLHS(self, GlobalLHS, LocalLHS, DoF):
		for row in range(24):
			for col in range(24):
				GlobalLHS[DoF[row], DoF[col]] += LocalLHS[row, col]
		return GlobalLHS

	def AssignDisplacementsToNodes(self, Nodes, DisplVector):  # Array of Nodes and the array of Displacements
		for node in range(self.__NumberOfNodes):
			Nodes[node].SetDisplacements(DisplVector[node*3], DisplVector[node*3 + 1], DisplVector[node*3 + 2])

	def Solve(Self, GlobalLHS, GlobalRHS):
		return np.linalg.solve(GlobalLHS, GlobalRHS)

	def CalculateAndAssignReactions(self, GlobalLHS, DisplVector, Nodes): # To check
		Reactions = np.dot(GlobalLHS, DisplVector)
		for node in range(self.__NumberOfNodes):
			Nodes[node].SetReactions(Reactions[node*3], Reactions[node*3 + 1], Reactions[node*3 + 2])
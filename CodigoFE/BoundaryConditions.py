import Node
import numpy as np

def Wait():
	input("Press Something")

class BoundaryCondition:
	def __init__(self, BoundaryConditionType, BoundaryContainer):  # BoundaryContainer = [Node, a,b,c, dx, dy, dz]
		self.__BoundaryConditionType = BoundaryConditionType  # String defining if its force or displ condition
		self.__BoundaryContainer  = BoundaryContainer
		self.__NumberOfFixedNodes = np.shape(BoundaryContainer)[0]

	def ApplyBoundaryCondition(self, GlobalLHS, GlobalRHS):
		if self.__BoundaryConditionType == "ImposedDisplacement":
			for node in range(self.__NumberOfFixedNodes):
				Node = self.__BoundaryContainer[node, 0]
				DoFNode = Node.GetNodeDoF()

				if self.__BoundaryContainer[node, 1] == 1:
					GlobalRHS[DoFNode[0]]    = self.__BoundaryContainer[node, 4]

					# Setting the row as identity matrix
					GlobalLHS[DoFNode[0], :] = np.zeros(np.shape(GlobalLHS)[0])
					GlobalLHS[DoFNode[0], DoFNode[0]] = 1.0

				if self.__BoundaryContainer[node, 2] == 1:
					GlobalRHS[DoFNode[1]] = self.__BoundaryContainer[node, 5]
					GlobalLHS[DoFNode[1], :] = np.zeros(np.shape(GlobalLHS)[0])
					GlobalLHS[DoFNode[1], DoFNode[1]] = 1.0					

				if self.__BoundaryContainer[node, 3] == 1:
					GlobalRHS[DoFNode[2]] = self.__BoundaryContainer[node, 6]
					GlobalLHS[DoFNode[2], :] = np.zeros(np.shape(GlobalLHS)[0])
					GlobalLHS[DoFNode[2], DoFNode[2]] = 1.0

		if self.__BoundaryConditionType == "NodalForce":
			for node in range(self.__NumberOfFixedNodes):
				Node = self.__BoundaryContainer[node, 0]
				DoFNode = Node.GetNodeDoF()

				if self.__BoundaryContainer[node, 1] == 1:
					GlobalRHS[DoFNode[0]] = self.__BoundaryContainer[node, 4]

				if self.__BoundaryContainer[node, 2] == 1:
					GlobalRHS[DoFNode[1]] = self.__BoundaryContainer[node, 5]				

				if self.__BoundaryContainer[node, 3] == 1:
					GlobalRHS[DoFNode[2]] = self.__BoundaryContainer[node, 6]
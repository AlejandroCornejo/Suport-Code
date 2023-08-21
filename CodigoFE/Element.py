import Node
import numpy as np

def Wait():
	input("Press Something")

class Element:
	def __init__(self, Id, Material, Node1, Node2, Node3, Node4, Node5, Node6, Node7, Node8):
		self.__Material = Material
		self.__Id = Id
		self.__Node1 = Node1
		self.__Node2 = Node2
		self.__Node3 = Node3
		self.__Node4 = Node4
		self.__Node5 = Node5
		self.__Node6 = Node6
		self.__Node7 = Node7
		self.__Node8 = Node8

		self.__NumberOfNodes = 8
		self.__StrainVectors = np.zeros((6, 8))  # Strain Vector in each Integration Point is a column vector
		self.__StressVectors = np.zeros((6, 8))
		self.__DeformationMatrices = np.zeros((6, 24, 8))

	def GetId(self):
		return self.__Id

	def GetNumberOfNodes(self):
		return self.__NumberOfNodes

	def GetNodes(self): # Returns an array of the nodes of thje element
		return [self.__Node1, self.__Node2, self.__Node3, self.__Node4, self.__Node5, self.__Node6, self.__Node7, self.__Node8]

	def SetDeformationMatrix(self, B, IntPoint):
		self.__DeformationMatrices[:,:,IntPoint] = B

	def GetDeformationMatrix(self, IntPoint):
		return self.__DeformationMatrices[:,:,IntPoint]

	def GetDoF(self):
		Nodes = self.GetNodes()
		DoF = [0.0]*24
		for node in range(self.GetNumberOfNodes()):
			DoFNode = Nodes[node].GetNodeDoF()
			DoF[node*3]     = DoFNode[0]
			DoF[node*3 + 1] = DoFNode[1]
			DoF[node*3 + 2] = DoFNode[2]
		return DoF

	def GetMaterial(self):
		return self.__Material

	def GetStrainVector(self, PointNumber):
		return self.__StrainVectors[:, PointNumber]

	def SetStrainVector(self, StrainVector, PointNumber):
		self.__StrainVectors[:, PointNumber] = StrainVector

	def GetStressVectors(self):
		return self.__StressVectors

	def GetStrainVectors(self):
		return self.__StrainVectors

	def GetStressVector(self, PointNumber):
		return self.__StressVectors[:, PointNumber]

	def SetStressVector(self, StressVector, PointNumber):
		self.__StressVectors[:, PointNumber] = StressVector

	def CalculateStrainVector(self, B):
		DisplVector = np.zeros(3 * self.GetNumberOfNodes())
		Nodes = self.GetNodes()

		for index in range(self.GetNumberOfNodes()):
			DisplVector[3 * index]     = Nodes[index].GetDisplacements()[0]
			DisplVector[3 * index + 1] = Nodes[index].GetDisplacements()[1]
			DisplVector[3 * index + 2] = Nodes[index].GetDisplacements()[2]

		return np.dot(B,DisplVector)

	def CalculateStressVector(self, StrainVector, C):
		return np.dot(C, StrainVector)

	def GetInfo(self):
		print("Linear Hexaedra Finite Element with 8 Nodes and 8 Integration Points")

	def GetLinearConstitutiveMatrix(self, E, nu):
		alpha = E*(1.0-nu)/(1.0+nu)/(1-2*nu)
		C = np.zeros((6,6))
		C = np.array([[1, nu/(1-nu), nu/(1-nu),0,0,0],
			[ nu/(1-nu),1, nu/(1-nu),0,0,0],
			[ nu/(1-nu), nu/(1-nu),1,0,0,0],
			[0,0,0,(1-2*nu)/(2*(1-nu)),0,0],
			[0,0,0,0,(1-2*nu)/(2*(1-nu)),0],
			[0,0,0,0,0,(1-2*nu)/(2*(1-nu))]])
		C *= alpha
		return C

	def GetLocalCoordinates(self):
		Node1 = Node.Node(1 ,-1, -1, -1)
		Node2 = Node.Node(2 , 1, -1, -1)
		Node3 = Node.Node(3 , 1,  1, -1)
		Node4 = Node.Node(4 ,-1,  1, -1)
		Node5 = Node.Node(5 ,-1, -1,  1)
		Node6 = Node.Node(6 , 1, -1,  1)
		Node7 = Node.Node(7 , 1,  1,  1)
		Node8 = Node.Node(8 ,-1,  1,  1)
		return [Node1, Node2, Node3, Node4, Node5, Node6, Node7, Node8]

	def GetLocalShapeFunctions(self, rPoint):  # Evaluates the Shape Functions at a given Point (x,y,z)
		ShapeValues = [0.0] * self.GetNumberOfNodes()

		ShapeValues[0] = 0.125*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 - rPoint[2] )
		ShapeValues[1] = 0.125*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 - rPoint[2] )
		ShapeValues[2] = 0.125*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 - rPoint[2] )
		ShapeValues[3] = 0.125*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 - rPoint[2] )
		ShapeValues[4] = 0.125*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 + rPoint[2] )
		ShapeValues[5] = 0.125*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 + rPoint[2] )
		ShapeValues[6] = 0.125*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 + rPoint[2] ) 
		ShapeValues[7] = 0.125*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 + rPoint[2] )
		return ShapeValues

	def GetLocalShapeFunctionsGradients(self, rPoint):  # [dNi/de_j ] so called DN_De
		rResult = np.zeros((self.GetNumberOfNodes() ,3))
		rResult[ 0, 0 ] = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] )  # DN1/Dchi
		rResult[ 0, 1 ] = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] )  # DN1/Deta
		rResult[ 0, 2 ] = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] )  # DN1/Dksi

		rResult[ 1, 0 ] =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] )
		rResult[ 1, 1 ] = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] )
		rResult[ 1, 2 ] = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] )

		rResult[ 2, 0 ] =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] )
		rResult[ 2, 1 ] =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] )
		rResult[ 2, 2 ] = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] )

		rResult[ 3, 0 ] = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] )
		rResult[ 3, 1 ] =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] )
		rResult[ 3, 2 ] = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] )

		rResult[ 4, 0 ] = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] )
		rResult[ 4, 1 ] = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] )
		rResult[ 4, 2 ] =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] )

		rResult[ 5, 0 ] =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] )
		rResult[ 5, 1 ] = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] )
		rResult[ 5, 2 ] =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] )

		rResult[ 6, 0 ] =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] )
		rResult[ 6, 1 ] =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] )
		rResult[ 6, 2 ] =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] )

		rResult[ 7, 0 ] = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] )
		rResult[ 7, 1 ] =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] )
		rResult[ 7, 2 ] =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] )

		return rResult

	def GetJacobianMatrix(self, IntegrationPoint): 
		Nodes = self.GetNodes() 
		JacobianMatrix = np.zeros((3,3))                
		CoordsIP  = IntegrationPoint.GetCoordinates()
		ShapeDeriv = self.GetLocalShapeFunctionsGradients(CoordsIP)  # Matrix 8 x 3
		
		for node in range(self.GetNumberOfNodes()):
			CoordNode = Nodes[node].GetCoordinates()

			JacobianMatrix[0,0] += ShapeDeriv[node,0]*CoordNode[0]
			JacobianMatrix[0,1] += ShapeDeriv[node,0]*CoordNode[1]
			JacobianMatrix[0,2] += ShapeDeriv[node,0]*CoordNode[2]

			JacobianMatrix[1,0] += ShapeDeriv[node,1]*CoordNode[0]
			JacobianMatrix[1,1] += ShapeDeriv[node,1]*CoordNode[1]
			JacobianMatrix[1,2] += ShapeDeriv[node,1]*CoordNode[2]

			JacobianMatrix[2,0] += ShapeDeriv[node,2]*CoordNode[0]
			JacobianMatrix[2,1] += ShapeDeriv[node,2]*CoordNode[1]
			JacobianMatrix[2,2] += ShapeDeriv[node,2]*CoordNode[2]

		return JacobianMatrix

	def CalculateDN_DX(self, DN_De, J): # Computes the cartesian derivatives of the Shape functions
		invJ = np.linalg.inv(J)
		DN_DX = np.zeros((self.GetNumberOfNodes() ,3))

		for node in range(self.GetNumberOfNodes()):
			M_DN_DX = np.dot(invJ, DN_De[node,:])
			DN_DX[node,0] = M_DN_DX[0]
			DN_DX[node,1] = M_DN_DX[1]
			DN_DX[node,2] = M_DN_DX[2]

		return DN_DX
		
	def CalculateDeformationMatrix(self, ShapeDeriv):   
		B = np.zeros((6,24))

		for node in range(self.GetNumberOfNodes()):
			B[0, node*3]     = ShapeDeriv[node,0]
			B[1, node*3 + 1] = ShapeDeriv[node,1]
			B[2, node*3 + 2] = ShapeDeriv[node,2]

			B[3, node*3]     = ShapeDeriv[node,1]
			B[3, node*3 + 1] = ShapeDeriv[node,0]

			B[4, node*3]     = ShapeDeriv[node,2]
			B[4, node*3 + 2] = ShapeDeriv[node,0]

			B[5, node*3 + 1] = ShapeDeriv[node,2]
			B[5, node*3 + 2] = ShapeDeriv[node,1]

		return B

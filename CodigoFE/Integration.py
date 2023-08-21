import Node
import scipy as sp 

class IntegrationQuadrature:
	def __init__(self):
		self.__IntegrationPoints = [0.0] * 8
		self.__IntegrationWeights = [0.0] * 8
		self.__NumberOfPoints = 8

	def SetIntegrationPoints(self):
			self.__IntegrationPoints[0] = Node.Node(1,-1.00/sp.sqrt(3.0) ,-1.00/sp.sqrt(3.0), -1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[1] = Node.Node(2 ,1.00/sp.sqrt(3.0), -1.00/sp.sqrt(3.0), -1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[2] = Node.Node(3 ,1.00/sp.sqrt(3.0) , 1.00/sp.sqrt(3.0), -1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[3] = Node.Node(4,-1.00/sp.sqrt(3.0) , 1.00/sp.sqrt(3.0), -1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[4] = Node.Node(5,-1.00/sp.sqrt(3.0) , -1.00/sp.sqrt(3.0), 1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[5] = Node.Node(6,1.00/sp.sqrt(3.0) ,  -1.00/sp.sqrt(3.0),  1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[6] = Node.Node(7, 1.00/sp.sqrt(3.0) ,  1.00/sp.sqrt(3.0), 1.00/sp.sqrt(3.0))
			self.__IntegrationPoints[7] = Node.Node(8,-1.00/sp.sqrt(3.0) ,  1.00/sp.sqrt(3.0), 1.00/sp.sqrt(3.0))

			for index in range(self.__NumberOfPoints):
				self.__IntegrationWeights[index] = 1.00

	def GetIntegrationPoints(self):
		points = [Node.Node(0,0,0,0)] * self.__NumberOfPoints
		for index in range(self.__NumberOfPoints):
			points[index] = self.__IntegrationPoints[index]
		return points

	def GetNumberOfIP(self):
		return self.__NumberOfPoints

	def GetIntegrationWeights(self):
		return self.__IntegrationWeights
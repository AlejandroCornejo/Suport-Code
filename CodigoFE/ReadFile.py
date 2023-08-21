import numpy as np
import Node
import Element
import Material

def Wait():
	input("Press Something")

class ReadFile:
	def __init__(self, InputFile):
		self.__LinesOfInput = InputFile.readlines()
		self.__NumberOfNodes      = int(self.__LinesOfInput[1].split("\t")[0])
		self.__NumberOfElements   = int(self.__LinesOfInput[1].split("\t")[1])
		self.__NumberOfMaterials  = int(self.__LinesOfInput[1].split("\t")[2])
		self.__ConnectivityMatrix = np.zeros((self.__NumberOfElements, 10))
		self.__NodesContainer     = np.zeros((self.__NumberOfNodes, 4))
		self.NodesWithBC         = 0
		self.NodesWithNodalForce = 0

	def GetNumberOfNodes(self):
		return self.__NumberOfNodes

	def GetNumberOfElements(self):
		return self.__NumberOfElements

	def GetNumberOfMaterials(self):
		return self.__NumberOfMaterials

	def GetMaterials(self):
		Materials = [0]*self.__NumberOfMaterials
		for nummat in range(self.__NumberOfMaterials):
			Id = int(self.__LinesOfInput[3+nummat].split("\t")[0])
			E  = float(self.__LinesOfInput[3+nummat].split("\t")[1])
			nu = float(self.__LinesOfInput[3+nummat].split("\t")[2])
			Materials[nummat] = Material.ElasticMaterial(Id, E, nu)
			
		return Materials

	def GetElements(self, Nodes, Materials):  # arrays of nodes and materials
		ElemCont = [0]*self.__NumberOfElements
		NumMat = self.__NumberOfMaterials
		NumElem = self.__NumberOfElements
		AuxIndex = int(0)

		for line in range(4 + NumMat, 4 + NumMat + NumElem):
			LineCont = self.__LinesOfInput[line].split("\t")
			LineCont = [int(i) for i in LineCont]

			ElemCont[AuxIndex] = Element.Element(LineCont[0], Materials[LineCont[1]-1], Nodes[LineCont[2]-1],
				Nodes[LineCont[3]-1],Nodes[LineCont[4]-1],Nodes[LineCont[5]-1],Nodes[LineCont[6]-1]
				,Nodes[LineCont[7]-1],Nodes[LineCont[8]-1],Nodes[LineCont[9]-1])
			AuxIndex += 1

		return ElemCont

	def GetNodes(self):
		NumElem = self.__NumberOfElements
		NumNodes = self.__NumberOfNodes
		NumMat = self.__NumberOfMaterials
		AuxIndex = int(0)
		Nodes = [0]*NumNodes

		for line in range(5+NumMat+NumElem, 5+NumMat+NumElem+NumNodes):
			LineCont = self.__LinesOfInput[line].split("\t")
			Nodes[AuxIndex] = Node.Node(int(LineCont[0]), float(LineCont[1]), float(LineCont[2]), float(LineCont[3]))
			AuxIndex += 1

		return Nodes

	def GetFixedDoF(self, Nodes):
		NumElem = self.__NumberOfElements
		NumNodes = self.__NumberOfNodes
		NumMat = self.__NumberOfMaterials
		LineCont = self.__LinesOfInput[6+NumMat+NumElem+NumNodes].split("\t")
		NumFixedNodes = int(LineCont[0])
		self.NodesWithBC = NumFixedNodes
		AuxIndex = 0

		FixedDoF = []
		for i in range(NumFixedNodes):
			FixedDoF.append([0]*7)

		for line in range(7+NumMat+NumElem+NumNodes, 7+NumMat+NumElem+NumNodes+NumFixedNodes):
			LineCont = self.__LinesOfInput[line].split("\t")

			FixedDoF[AuxIndex][0] = Nodes[int(LineCont[0])-1]
			FixedDoF[AuxIndex][1] = int(LineCont[1])
			FixedDoF[AuxIndex][2] = int(LineCont[2])
			FixedDoF[AuxIndex][3] = int(LineCont[3])
			FixedDoF[AuxIndex][4] = float(LineCont[4])
			FixedDoF[AuxIndex][5] = float(LineCont[5])
			FixedDoF[AuxIndex][6] = float(LineCont[6])
			AuxIndex += 1

		FixedDoF  = np.array(FixedDoF)
		return FixedDoF

	def GetNodalForces(self, Nodes):
		NumElem = self.__NumberOfElements
		NumNodes = self.__NumberOfNodes
		NumMat = self.__NumberOfMaterials
		NumFixedNodes = self.NodesWithBC
		LineCont = self.__LinesOfInput[8+NumMat+NumElem+NumNodes+NumFixedNodes].split("\t")
		NumNodalForces = int(LineCont[0])
		AuxIndex = 0

		FixedDoF = []
		for i in range(NumFixedNodes):
			FixedDoF.append([0]*7)

		for line in range(9+NumMat+NumElem+NumNodes+NumFixedNodes, 9+NumMat+NumElem+NumNodes+NumFixedNodes+NumNodalForces):
			LineCont = self.__LinesOfInput[line].split("\t")

			FixedDoF[AuxIndex][0] = Nodes[int(LineCont[0])-1]
			FixedDoF[AuxIndex][1] = int(LineCont[1])
			FixedDoF[AuxIndex][2] = int(LineCont[2])
			FixedDoF[AuxIndex][3] = int(LineCont[3])
			FixedDoF[AuxIndex][4] = float(LineCont[4])
			FixedDoF[AuxIndex][5] = float(LineCont[5])
			FixedDoF[AuxIndex][6] = float(LineCont[6])
			AuxIndex += 1

		FixedDoF  = np.array(FixedDoF)
		return FixedDoF	
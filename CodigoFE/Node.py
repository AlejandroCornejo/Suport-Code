import scipy as sp
import numpy as np

class Node:
	def __init__(self, Id, Xcoord, Ycoord, Zcoord):
		self.__Id = Id
		self.__Xcoord = Xcoord
		self.__Ycoord = Ycoord
		self.__Zcoord = Zcoord

		self.__Xdispl = 0.0
		self.__Ydispl = 0.0
		self.__Zdispl = 0.0

		self.__XReac = 0.0
		self.__YReac = 0.0
		self.__ZReac = 0.0

	def GetId(self):
		return self.__Id

	def GetX(self):
		return self.__Xcoord

	def GetY(self):
		return self.__Ycoord		

	def GetZ(self):
		return self.__Zcoord

	def GetNodeDoF(self):
		return [(self.GetId()-1)*3, (self.GetId()-1)*3 + 1, (self.GetId()-1)*3 + 2]

	def GetCoordinates(self):
		return [self.__Xcoord, self.__Ycoord, self.__Zcoord]

	def SetDisplacements(self, Xdispl, Ydispl, Zdispl):
		self.__Xdispl = Xdispl
		self.__Ydispl = Ydispl
		self.__Zdispl = Zdispl

	def SetReactions(self, Xdispl, Ydispl, Zdispl):
		self.__XReac = Xdispl
		self.__YReac = Ydispl
		self.__ZReac = Zdispl

	def GetDisplacements(self):
		return [self.__Xdispl, self.__Ydispl, self.__Zdispl]

	def GetReactions(self):
		return [self.__XReac, self.__YReac, self.__ZReac]

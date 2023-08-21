
class ElasticMaterial:
	def __init__(self, Id = 1, Young = 35000 , Poisson = 0.22):
		self.__Id = Id
		self.__Young = Young 
		self.__Poisson = Poisson 
		self.__ConstitutiveLaw = "LinearElastic"

	def GetYoungModulus(self):
		return self.__Young

	def GetPoisson(self):
		return self.__Poisson

from . import grayencoder as ge

class Phenome:
	def __init__(self,ind):
		self.ind = ind
		self.genome = self.constructGenome(self.ind)
		self.constructPhenome(self.genome)

	def constructGenome(self,ind):
		#genome construction goes here for inherited classes
		return None

	def constructPhenome(self,genome):
		#phenome construction goes here for inherited classes
		return None


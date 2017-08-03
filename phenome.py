import grayencoder as ge

class Phenome:
	def __init__(self,ind):
		self.ind = ind
		self.genome = self.constructGenome(self.ind)
		self.constructPhenome(self.genome)

	def constructGenome(self,ind):
		test = ge.genCode(5)
		print(test)
		print(ge.readCode(test))

	def constructPhenome(self,genome):
		return


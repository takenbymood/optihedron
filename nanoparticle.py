from GA import phenome
from GA import grayencoder as ge
from Tools import listtools as lt

class Ligand:
	def __init__(self,eps,sig,rad,ang,mass=1,cutoff=2.5):
		self.rad = rad
		self.ang = ang
		self.eps = eps
		self.sig = sig
		self.mass = mass
		self.cutoff = cutoff

	def __str__(self):
		ligstr= "rad:"+str(self.rad)+", ang:"+str(self.ang)+", eps:"+str(self.eps)+", sig:"+str(self.sig)
		ligstr+= ", mass:"+str(self.mass)+", cutoff:"+str(self.cutoff)
		return ligstr

class NanoParticle:
	ligands = []
	def __init__(
		self,
		x=0,
		y=0,
		mass=1,
		eps=1,
		sig=4,
		cutoff=2**(1/6)
		):
		self.x = x
		self.y = y
		self.mass = mass
		self.eps = eps
		self.sig = sig
		self.cutoff = cutoff
		self.ligands=[]

	def addLigand(self,ligand):
		if(isinstance(ligand,Ligand)):
			self.ligands.append(ligand)

	def spaceIsOccupied(self,targetAngle):
		for l in self.ligands:
			if targetAngle < l.ang + 0.7 and targetAngle > l.ang - 0.7:
				return True
		return False

	def findNearestSpace(self,targetAngle,spacing):
		freeSpace = targetAngle
		if self.spaceIsOccupied(targetAngle):
			return(self.findNearestSpace(targetAngle+float(random.randint(0,2)-1)*spacing,spacing))
		else:
			return freeSpace

	def __str__(self):
		protstr = "x:"+str(self.x)+", y:"+str(self.y)+", m:"+str(self.mass)
		protstr += ", eps:"+str(self.eps)+", sig:"+str(self.sig)+", cut:"+str(self.cutoff)
		protstr += "\n"+str(len(self.ligands))+" ligands"
		i=0
		for l in self.ligands:
			i+=1
			protstr += "\nligand " + str(i) +" - "+ str(l)
		return protstr

class NanoParticlePhenome(phenome.Phenome):

	def __init__(
		self,
		ind,
		epsPlaces,
		angPlaces,
		minEps,
		maxEps
		):
		self.epsPlaces = epsPlaces
		self.angPlaces = angPlaces
		self.minEps = minEps
		self.maxEps = maxEps
		phenome.Phenome.__init__(self,ind)

	def constructGenome(self,ind):
		geneSize = self.epsPlaces+self.angPlaces
		if len(ind) < geneSize:
			print("genome too short to construct phenome")
			return None
		genelist = lt.subdivide(list(ind),geneSize)
		genes = []
		for g in genelist:
			gene = {}
			epsGene = g[0:self.epsPlaces]
			angGene = g[self.epsPlaces:self.epsPlaces+self.angPlaces]
			gene['eps'] = (ge.read(epsGene)*(self.maxEps-self.minEps))/ge.max(epsGene)
			gene['ang'] = (ge.read(angGene)*(6.2831)/ge.max(angGene))
			genes.append(gene)
		return genes


	def constructPhenome(self,ind):
		self.particle = NanoParticle()
		for g in self.genome:
			if not self.particle.spaceIsOccupied(g['ang']):
				self.particle.addLigand(Ligand(g['eps'],1,4,g['ang']))
		return self.particle
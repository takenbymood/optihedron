from ga import phenome
from ga import grayencoder as ge
from tools import listtools as lt
from tools.vectools import polarToCartesianVector 
from tools import icosatiler
import numpy as np

class Ligand:
	def __init__(self,eps,sig,rad,polAng,aziAng,mass=1,cutoff=1.8):
		self.rad = rad
		self.polAng = polAng
		self.aziAng = aziAng
		self.eps = eps
		self.sig = sig
		self.mass = mass
		self.cutoff = cutoff
		self.size = 1.0

	def __str__(self):
		ligstr= "rad:"+str(self.rad)+", polAng:"+str(self.polAng)+", aziAng:"+str(self.aziAng)+", eps:"+str(self.eps)
		ligstr+= ", sig:"+str(self.sig)+", mass:"+str(self.mass)+", cutoff:"+str(self.cutoff)
		return ligstr

class NanoParticle:
	ligands = []
	def __init__(
		self,
		mass=1.0,
		eps=1.0,
		sig=4.0,
		cutoff=2.0**(1.0/6.0)
		):
		self.mass = mass
		self.eps = eps
		self.sig = sig
		self.cutoff = cutoff
		self.ligands=[]

	def addLigand(self,ligand):
		if(isinstance(ligand,Ligand)):
			self.ligands.append(ligand)

	def spaceIsOccupied(self,polarTargetAngle,azimuthalTargetAngle):
		target_v = polarToCartesianVector(self.sig, polarTargetAngle, azimuthalTargetAngle)
		for ligand in self.ligands:
			ligand_v = polarToCartesianVector(ligand.rad, ligand.polAng, ligand.aziAng)
			d = np.linalg.norm(np.subtract(ligand_v,target_v))
			if d<ligand.size:
				return True
		return False

	def __str__(self):
		protstr = "m:"+str(self.mass)+", eps:"+str(self.eps)+", sig:"+str(self.sig)+", cut:"+str(self.cutoff)
		protstr += "\n"+str(len(self.ligands))+" ligands"
		for i,l in zip(range(len(self.ligands)),self.ligands):
			protstr += "\nligand " + str(i+1) +" - "+ str(l)
		return protstr

class NanoParticlePhenome(phenome.Phenome):

	def __init__(
		self,
		ind,
		exprPlaces,
		epsPlaces,
		polarAngPlaces,
		azimuthalAngPlaces,
		minEps,
		maxEps
		):
		self.exprPlaces = exprPlaces
		self.epsPlaces = epsPlaces
		self.polarAngPlaces = polarAngPlaces
		self.azimuthalAngPlaces = azimuthalAngPlaces
		self.minEps = minEps
		self.maxEps = maxEps
		self.radius = 4
		phenome.Phenome.__init__(self,ind)

	def constructGenelist(self,ind):
		geneSize = self.exprPlaces+self.epsPlaces+self.polarAngPlaces+self.azimuthalAngPlaces
		if len(ind) < geneSize:
			print("genome too short to construct phenome")
			return None
		genelist = lt.subdivide(list(ind),geneSize)
		if len(genelist[-1]) < geneSize:
			genelist = genelist[:-1]
		return genelist


	def constructGenome(self,ind):
		genes = []
		for g in self.genelist:
			gene = {}
			exprGene = g[0:self.exprPlaces]			
			epsGene = g[self.exprPlaces:self.exprPlaces+self.epsPlaces]
			polarAngGene = g[self.exprPlaces+self.epsPlaces:self.exprPlaces+self.epsPlaces+self.polarAngPlaces]
			azimuthalAngGene = g[self.exprPlaces+self.epsPlaces+self.polarAngPlaces:self.exprPlaces+self.epsPlaces+self.polarAngPlaces+self.azimuthalAngPlaces]
			gene['expr'] = True if ge.read(exprGene)/ge.max(exprGene) >= 0.5 else False 			
			gene['eps'] = ((ge.read(epsGene)*(self.maxEps-self.minEps))/ge.max(epsGene))+self.minEps
			gene['polAng'] = (ge.read(polarAngGene)*(3.141)/ge.max(polarAngGene))
			gene['aziAng'] = (ge.read(azimuthalAngGene)*(6.2831)/ge.max(azimuthalAngGene))
			genes.append(gene)
		return genes


	def constructPhenome(self,ind):
		self.particle = NanoParticle()
		for g in self.genome:
			if g['expr']:				
				if not self.particle.spaceIsOccupied(g['polAng'],g['aziAng']):
					self.particle.addLigand(Ligand(g['eps'],1,self.radius,g['polAng'],g['aziAng']))			
		return self.particle

class CoveredNanoParticlePhenome(phenome.Phenome):
	def __init__(
		self,
		ind,
		exprPlaces,
		epsPlaces,
		minEps,
		maxEps
		):
		self.exprPlaces = exprPlaces
		self.epsPlaces = epsPlaces
		self.minEps = minEps
		self.maxEps = maxEps
		self.radius = 4
		phenome.Phenome.__init__(self,ind)

	def constructGenelist(self,ind):
		geneSize = self.exprPlaces+self.epsPlaces
		if len(ind) < geneSize:
			print("genome too short to construct phenome")
			return None
		genelist = lt.subdivide(list(ind),geneSize)
		if len(genelist[-1]) < geneSize:
			genelist = genelist[:-1]
		return genelist

	def constructGenome(self,ind):
		genes = []
		for g in self.genelist:
			gene = {}
			exprGene = g[0:self.exprPlaces]			
			epsGene = g[self.exprPlaces:self.exprPlaces+self.epsPlaces]
			gene['expr'] = True if ge.read(exprGene)/ge.max(exprGene) >= 0.5 else False 			
			gene['eps'] = ((ge.read(epsGene)*(self.maxEps-self.minEps))/ge.max(epsGene))+self.minEps
			genes.append(gene)
		return genes

	def constructPhenome(self,ind):
		ligandPositions = icosatiler.cover72SpherePolar(self.radius)
		self.particle = NanoParticle()
		i = 0
		genomeOverrun = False;
		for g in self.genome:
			if i<len(ligandPositions):	
				v = ligandPositions[i]		
				eps = g['eps'] if g['expr'] else 0
				self.particle.addLigand(Ligand(eps,1,v[0],v[1],v[2]))
			elif i>=len(ligandPositions):
				genomeOverrun = True
			i+=1	
		if genomeOverrun:
			print('genome was too long for the symmetry used')		
		return self.particle

	


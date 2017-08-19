import sys
import os
import math

from Sims import lammpsbuilder as lb

class MembraneSimulation(lb.LammpsSimulation):

	def __init__(
		self,
		name,
		protein,
		run,
		timestep,
		outdir,
		filedir,
		membraneFile,
		corepos_x=0, 
		corepos_y=10,
		dumpres="100"
		):
		lb.LammpsSimulation.__init__(
			self,
			name,
			run,
			filedir
			)
		self.setScript(lb.LammpsScript(
			run,
			timestep,
			"2",
			"xyz",
			"100",
			"out.xyz",
			'300',
			'0',
			'all create 1.0 1000',
			'lj/cut 5.04',
			'molecular',
			read_data=self.dataName
			))
		self.setData(lb.LammpsData(
			-75,
			75,
			-50,
			50,
			-1,
			1
			))
		self.outdir = outdir
		self.script.dump = "id all xyz "+dumpres+" "+outdir+name +"_out.xyz"
		self.data.atomTypes = 3+len(protein.ligands)
		self.data.bondTypes = 1
		self.data.angleTypes = 1
		self.data.addMass(1,1)
		self.data.addMass(2,1)
		self.data.addMass(3,1)
		for i in range(len(protein.ligands)):
			self.data.addMass(4+i,1)

		#startX = -(0.5*mLength*spacing)

		#self.data.addAtom(2,startX,0)

		# for i in range(mLength-2):
		# 	self.data.addAtom(1,startX+spacing*i+2,0)
		# 	self.data.addBond(1,i+1,i+2)
		# 	self.data.addAngle(1,i+1,i+2,i+3)
		if os.path.exists(membraneFile):
			a = self.data.addXyzFile(membraneFile)

		#self.data.addAtom(2,startX+spacing*mLength,0)
		#self.data.addBond(1,mLength-1,mLength)

		for i in range(2,a[1]+1):
			pAtom = a[0]+i
			self.data.addBond(1,pAtom-1,pAtom)

		self.data.addBond(1,1,a[1])

		for i in range(3,a[1]+1):
			pAtom = a[0]+i
			self.data.addAngle(1,pAtom-2,pAtom-1,pAtom)

		self.data.addAngle(1,a[1]-1,a[1],1)
		self.data.addAngle(1,a[1],1,2)

		mol = self.data.addAtom(3,corepos_x,corepos_y,0)

		self.script.addBond(1,30,1.3)
		self.script.addAngle(1,30,180)
		self.script.addPair("*","*",0,0,0)

		aType = 4
		for l in protein.ligands:
			self.data.addAtom(aType,corepos_x+l.rad*math.cos(l.ang),corepos_y+l.rad*math.sin(l.ang),0,mol)
			self.script.addPair("1",str(aType),l.eps,l.sig,l.sig*l.cutoff)
			aType+=1
		
		self.script.addPair(1,3,100,4,4.45)
		self.script.addPair(1,1,100,1.0,1.12)
		self.script.addPair(1,2,100,1.0,1.12)
		self.script.addPairModify("shift yes")

		self.script.addGroup("move",[1])
		self.script.addGroup("anchor",[2])
		self.script.addGroup("lipid",[1,2])
		pGroup = [3]
		for i in range(len(protein.ligands)):
			pGroup.append(i+4)
		self.script.addGroup("protein",pGroup)
		self.script.addPreFixLine("velocity move create 1.0 1")
		self.script.addPreFixLine("velocity protein set 0 -2 0")
		self.script.addFix("all","enforce2d")
		self.script.addFix("lipid","nph x 0.0 0.0 1.0 y 0.0 0.0 1.0 couple xy")
		self.script.addFix("protein","rigid/nve molecule")
		self.script.addFix("all","langevin 1 1 1 1000")
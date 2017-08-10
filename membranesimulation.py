class MembraneSimulation(lb.LammpsSimulation):

	def __init__(self,name,protein,outdir,filedir,mLength=120,spacing=1.3,corepos_x=0, corepos_y=10,run="200000",dumpres="100"):
		lb.LammpsSimulation.__init__(self,name,"out/",run=run)
		self.script.dump = "id all xyz "+dumpres+" "+outdir+name +"_out.xyz"
		self.data.atomTypes = 3+len(protein.ligands)
		self.data.bondTypes = 1
		self.data.angleTypes = 1
		self.data.addMass(1,1)
		self.data.addMass(2,1)
		self.data.addMass(3,3)
		for i in range(len(protein.ligands)):
			self.data.addMass(4+i,0.01)

		startX = -(0.5*mLength*spacing)

		self.data.addAtom(2,startX,0)

		for i in range(mLength-2):
			self.data.addAtom(1,startX+spacing*i+2,0)
			self.data.addBond(1,i+1,i+2)
			self.data.addAngle(1,i+1,i+2,i+3)

		self.data.addAtom(2,startX+spacing*mLength,0)
		self.data.addBond(1,mLength-1,mLength)

		mol = self.data.addAtom(3,corepos_x,corepos_y,0)

		self.script.addBond(1,2.0,1.3)
		self.script.addAngle(1,30,180)
		self.script.addPair("*","*",0,0,0)

		aType = 4
		for l in protein.ligands:
			self.data.addAtom(aType,corepos_x+l.rad*math.cos(l.ang),corepos_y+l.rad*math.sin(l.ang),0,mol)
			self.script.addPair("1",str(aType),l.eps,l.sig,l.sig*l.cutoff)
			aType+=1
		
		self.script.addPair(1,3,100,4.5,5.0)
		self.script.addPair(1,1,100,1.0,1.12246)
		self.script.addPair(1,2,100,1.0,1.12246)
		self.script.addPairModify("shift yes")

		self.script.addGroup("move",[1])
		self.script.addGroup("anchor",[2])
		pGroup = [3]
		for i in range(len(protein.ligands)):
			pGroup.append(i+4)
		self.script.addGroup("protein",pGroup)
		self.script.addFix("all","enforce2d")
		self.script.addFix("move","nve")
		self.script.addFix("protein","rigid/nve molecule")
		self.script.addFix("all","langevin 1 1 1 1000")
		
		self.script.addLine("fix 4 all wall/lj93 yhi 18 1.0 1.0 1.12")
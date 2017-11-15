
import sys
import os


class LammpsAtom:

	def __init__(
		self,
		atomId,
		moleculeId,
		atomType,
		x=0,
		y=0,
		z=0
		):
		self.atomId = atomId
		self.moleculeId = moleculeId
		self.x = x
		self.y = y
		self.z = z
		self.atomType = atomType

	def __str__(self):
		s=str(self.atomId)+" "+str(self.moleculeId)+" "+str(self.atomType)
		s+=" "+str(self.x)+" "+str(self.y)+" "+str(self.z)
		return s

class LammpsBond:

	def __init__(self,bondId,bondType,atom1,atom2):
		self.bondId = bondId
		self.bondType = bondType
		self.atom1 = atom1
		self.atom2 = atom2

	def __str__(self):
		s=str(self.bondId)+" "+str(self.bondType)+" "+str(self.atom1)+" "+str(self.atom2)
		return(s)

class LammpsAngle:

	def __init__(self,angleId,angleType,atom1,atom2,atom3):
		self.angleId = angleId
		self.angleType = angleType
		self.atom1 = atom1
		self.atom2 = atom2
		self.atom3 = atom3

	def __str__(self):
		s = str(self.angleId)+" "+str(self.angleType)+" "+str(self.atom1)+" "+str(self.atom2)+" "+str(self.atom3)
		return s

class LammpsMass:
	def __init__(self,atomType,mass):
		self.atomType = atomType
		self.mass = mass

	def __str__(self):
		s = str(self.atomType)+" "+str(self.mass)
		return s

class LammpsData:


	def __init__(self,
		xlo,
		xhi,
		ylo,
		yhi,
		zlo,
		zhi,
		atomTypes=0,
		bondTypes=0,
		angleTypes=0,
		):
		self.atomTypes = atomTypes
		self.bondTypes = bondTypes
		self.angleTypes = angleTypes
		self.xlo = xlo
		self.xhi = xhi
		self.ylo = ylo
		self.yhi = yhi
		self.zlo = zlo
		self.zhi = zhi
		self.atoms = []
		self.bonds = []
		self.masses = []
		self.angles = []

	def addXyzFile(self,filename):
		nAtoms = 0
		startId = len(self.atoms)
		if os.path.exists(filename):
			with open(filename, 'rb') as f:
				try:
					content = f.readlines()
					for l in content:
						atomData = l.replace('\n','').split(' ')
						nAtoms+=1
						self.addAtom(int(atomData[0]),float(atomData[1]),float(atomData[2]),float(atomData[3]))
				except : 
					print(str(filename) + " does not exist")
		return (startId,nAtoms)

	def addAtom(self,atomType,x,y,z=0,moleculeId=-1):
		atomId = len(self.atoms)+1
		if(moleculeId == -1):
			moleculeId = atomId
		a = LammpsAtom(atomId,moleculeId,atomType,x,y,z)
		self.atoms.append(a)
		return atomId

	def addBond(self,bondType,atom1,atom2):
		bondId = len(self.bonds)+1
		bond = LammpsBond(bondId,bondType,atom1,atom2)
		self.bonds.append(bond)
		return bondId

	def addAngle(self,angleType,atom1,atom2,atom3):
		angleId = len(self.angles)+1
		angle = LammpsAngle(angleId,angleType,atom1,atom2,atom3)
		self.angles.append(angle)
		return angleId

	def addMass(self,atomType,mass):
		self.masses.append(LammpsMass(atomType,mass))

	def buildFile(self):
		s = "LAMMPS CONFIG FILE\n\n"

		s += "  "+str(len(self.atoms))+" atoms\n"
		s += "  "+str(len(self.bonds))+" bonds\n"
		s += "  "+str(len(self.angles))+" angles\n"
		s += "\n"
		s += "  "+str(self.atomTypes)+" atom types\n"
		s += "  "+str(self.bondTypes)+" bond types\n"
		s += "  "+str(self.angleTypes)+" angle types\n"
		s += "\n"
		s += "  "+str(self.xlo)+"  "+str(self.xhi)+" xlo xhi\n"
		s += "  "+str(self.ylo)+"  "+str(self.yhi)+" ylo yhi\n"
		s += "  "+str(self.zlo)+"  "+str(self.zhi)+" zlo zhi\n"
		s += "\n"
		if(len(self.masses)>0):
			s+="Masses\n\n"
			for m in self.masses:
				s+= "  "+str(m)+"\n"
			s+="\n"
		if(len(self.atoms)>0):
			s+="Atoms\n\n"
			for a in self.atoms:
				s+="  "+str(a)+"\n"
			s+="\n"
		if(len(self.bonds)>0):
			s+="Bonds\n\n"
			for b in self.bonds:
				s+="  "+str(b)+"\n"
			s+="\n"
		if(len(self.angles)>0):
			s+="Angles\n\n"
			for a in self.angles:
				s+="  "+str(a)+"\n"
			s+="\n"
		return s

	def __str__(self):
		return self.buildFile()


class LammpsScript:
	sep = "			"

	def addPair(self,atom1,atom2,eps=0,sig=0,cutoff=""):
		s = str(atom1)+" "+str(atom2)+" "+str(eps)+" "+str(sig)+" "+str(cutoff)
		self.pair_coeffs.append(s)

	def addPairModify(self,line):
		self.pairmod.append(line)

	def addBond(self,bond,K,x0):
		s = str(bond)+" "+str(K)+" "+str(x0)
		self.bond_coeffs.append(s)

	def addAngle(self,angle,K,theta0):
		s = str(angle)+" "+str(K)+" "+str(theta0)
		self.angle_coeffs.append(s)

	def addGroup(self,name,members,order="type"):
		s = str(name)+" "+str(order)+" "
		for m in members:
			s+=str(m)+" "
		self.groups.append(s)

	def addFix(self,group,action):
		s = str(group)+ " " + str(action)
		self.fixes.append(s)

	def addPostFixLine(self,line):
		self.postLines.append(str(line))

	def addPreFixLine(self,line):
		self.preLines.append(str(line))

	def __init__(
		self,
		run,
		timestep,
		dimension,
		dumpStyle,
		dumpStep,
		dumpFile,
		thermo,
		seed,
		velocity,
		pair_style,
		atom_style,
		angle_style="harmonic",
		bond_style="harmonic",
		atom_modify="sort 0 1",
		neighbour="0.3 bin",
		neigh_modify="every 1 delay 1",
		units="lj",
		read_data=""
		):
		self.read_data = read_data
		self.dump = "id all " + dumpStyle + " " + dumpStep + " " + dumpFile
		self.dimension = str(dimension)
		self.units = units
		self.atom_style = atom_style
		self.atom_modify = atom_modify
		self.neighbour = neighbour
		self.neigh_modify = neigh_modify
		self.angle_style = angle_style
		self.bond_style = bond_style
		self.pair_style = pair_style
		self.timestep = str(timestep)
		self.thermo = str(thermo)
		self.run = str(run)
		self.velocity = "all create " + velocity + " " + str(seed)
		self.bond_coeffs=[]
		self.pair_coeffs=[]
		self.angle_coeffs=[]
		self.groups = []
		self.fixes = []
		self.preLines = []
		self.postLines = []
		self.pairmod=[]

	def validate(self):
		#check for missing properties
		return True

	def buildHeader(self):
		d = self.sep
		s="dimension"+d+self.dimension+"\n"
		s+="units"+d+self.units+"\n"
		s+="atom_style"+d+self.atom_style+"\n"
		s+="boundary p p p\n"
		s+="atom_modify"+d+self.atom_modify+"\n"
		s+="\n"
		s+="read_data"+d+self.read_data+"\n"
		#s+="neighbour"+d+self.neighbour+"\n"
		s+="neigh_modify"+d+self.neigh_modify+"\n"
		return s

	def buildCoefficients(self):
		d = self.sep
		s = ""
		if len(self.bond_coeffs)>0:
			s+="bond_style"+d+self.bond_style+"\n"
		for b in self.bond_coeffs:
			s+="bond_coeff"+d+b+"\n"
		if len(self.angle_coeffs)>0:
			s+="angle_style"+d+self.angle_style+"\n"
		for a in self.angle_coeffs:
			s+="angle_coeff"+d+a+"\n"
		if len(self.pair_coeffs)>0:
			s+="pair_style"+d+self.pair_style+"\n"
		for p in self.pair_coeffs:
			s+="pair_coeff"+d+p+"\n"
		for p in self.pairmod:
			s+="pair_modify"+d+p+"\n"
		return s

	def buildDump(self):
		d = self.sep
		s = ""
		s+="dump"+d+self.dump+"\n"
		return s

	def buildGroups(self):
		d = self.sep
		s = ""
		for g in self.groups:
			s+="group"+d+g+"\n"
		return s

	def buildPreFixLines(self):
		d = self.sep
		s=""
		for l in self.preLines:
			s+=l
			s+="\n"
		return s

	def buildFixes(self):
		d = self.sep
		s=""
		i=0
		for f in self.fixes:
			s+="fix"+d+str(i)+" "+f+"\n"
			i+=1
		s+="\n"
		return s

	def buildPostFixLines(self):
		d = self.sep
		s=""
		for l in self.postLines:
			s+=l
			s+="\n"
		return s

	def buildFooter(self):
		d = self.sep
		s=""
		s+="thermo"+d+self.thermo+"\n"
		s+="timestep"+d+self.timestep+"\n"
		s+="run"+d+str(self.run)+"\n"
		return s

	def buildFile(self):
		d = self.sep
		s = self.buildHeader()
		s += self.buildCoefficients()
		s += self.buildDump()
		s += self.buildGroups()
		s += self.buildPreFixLines()
		s += self.buildFixes()
		s += self.buildPostFixLines()
		s += self.buildFooter()
		return s	

	def __str__(self):
		return self.buildFile()
		

class LammpsSimulation:

	def __init__(self,name,run,filedir=""):
		self.name = name
		self.scriptName = name+"_script.in"
		self.dataName = name+"_data.data"
		self.filedir = filedir
		self.script = None
		self.data = None

	def __del__(self):
		name = ""
		scriptName = ""
		dataName = ""
		filedir = ""
		script = None
		data = None

	def setScript(self,script):
		self.script = script

	def setData(self,data):
		self.data = data

	def saveFiles(self):
		if self.script != None:
			with open(os.path.join(self.filedir,self.scriptName), 'w') as file_:
				file_.write(str(self.script))

		if self.data != None:
			with open(os.path.join(self.filedir,self.dataName), 'w') as file_:
				file_.write(str(self.data))

	def deleteFiles(self):
		os.remove(os.path.join(self.filedir,self.scriptName))
		os.remove(os.path.join(self.filedir,self.dataName))

	def __str__(self):
		return self.name + "\n" + str(self.script) + "\n" + str(self.data)

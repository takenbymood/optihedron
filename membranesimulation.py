import os
import shutil
import math
import numpy as np
from tools import templatetools as tt
from tools import vectools
from tools import misctools
import random

class MembraneSimulation():
	def __init__(
		self,
		name, 
		protein, 
		run, 
		timestep, 
		outdir, 
		filedir, 
		datatemplate, 
		scripttemplate, 
		corepos_x=0, 
		corepos_y=0, 
		corepos_z=6.5, 
		dumpres="100",
		rAxis = [0,0,1],
		rAmount = 0.0
		):

		self.name = name
		self.scriptName = name+'_script.in'
		self.dataName = name+"_data.data"
		self.outName = name+"_out.xyz"
		self.protein = protein
		self.run = run
		self.timestep = timestep		
		self.outdir = outdir
		self.filedir = filedir
		self.datatemplate = datatemplate
		self.scripttemplate = scripttemplate
		self.corepos_x = corepos_x
		self.corepos_y = corepos_y
		self.corepos_z = corepos_z
		self.dumpres = dumpres
		self.nonLigandAtomCount = 2900 + 1
		self.rAxis = rAxis
		self.rAmount = rAmount
		self.rmat = vectools.buildERMatrix(self.rAxis, self.rAmount)

	def saveFiles(self):
		scratch = os.path.join(self.filedir, self.name)		

		simData = os.path.join(self.filedir, self.dataName)
		simScript = os.path.join(self.filedir, self.scriptName)
		shutil.copyfile(self.datatemplate, simData)
		shutil.copyfile(self.scripttemplate, simScript)

		ligandMasses = ''		
		npPositions = '{0} 2 {1} {2} {3}   1 1 0   0 0 0\n'.format(self.nonLigandAtomCount, self.corepos_x, self.corepos_y, self.corepos_z)
		npVelocities = '{0} 0 0 0 0 0 0\n'.format(self.nonLigandAtomCount)
		ligandMembraneInteractions = ''
		

		for i, ligand in enumerate(self.protein.ligands, 1):
			#physics convention, theta = polar, phi = azimuthal
			ligand_v = vectools.polarToCartesianVector(ligand.rad, ligand.polAng, ligand.aziAng)

			r_ligand_v = np.dot(self.rmat,ligand_v)

			ligand_x = self.corepos_x+r_ligand_v[0]
			ligand_y = self.corepos_y+r_ligand_v[1]
			ligand_z = self.corepos_z+r_ligand_v[2]

			ligandMasses += '{0} {1}\n'.format(2+i, ligand.mass)		
			npPositions += '{0} {1} {2} {3} {4}  1 1 0   0 0 0\n'.format(self.nonLigandAtomCount+i,2+i,ligand_x,ligand_y,ligand_z)
			npVelocities += '{0} 0 0 0 0 0 0\n'.format(self.nonLigandAtomCount+i)
			ligandMembraneInteractions += 'pair_coeff		1	{}	lj/cut		{}	{}	{}\n'.format(2+i,ligand.eps,ligand.sig,ligand.sig*ligand.cutoff)

		dataTemp = tt.loadTemplate(self.datatemplate)
		scriptTemp = tt.loadTemplate(self.scripttemplate)
		dataTemp = dataTemp.replace('_ATOM COUNT PLACEHOLDER_', '{} atoms\n'.format(self.nonLigandAtomCount+len(self.protein.ligands)))
		dataTemp = dataTemp.replace('_ATOM TYPE COUNT PLACEHOLDER_', '{} atom types\n'.format(2+len(self.protein.ligands)))
		dataTemp = dataTemp.replace('_LIGAND MASSES PLACEHOLDER_', ligandMasses)
		dataTemp = dataTemp.replace('_NANOPARTICLE POSITIONS PLACEHOLDER_', npPositions)
		dataTemp = dataTemp.replace('_NANOPARTICLE VELOCITIES PLACEHOLDER_', npVelocities)

		scriptTemp = scriptTemp.replace('_DATA FILE PLACEHOLDER_', 'read_data			"{}"\n'.format(simData))
		scriptTemp = scriptTemp.replace('_LIGAND GROUP PLACEHOLDER_', 'group				ligand 	type {}:{}\n'.format(3,2+len(self.protein.ligands)))
		scriptTemp = scriptTemp.replace('_NANOPARTICLE GROUP PLACEHOLDER_', 'group				np      type 2:{}\n'.format(2+len(self.protein.ligands)))
		scriptTemp = scriptTemp.replace('_LIGAND MEMBRANE INTERACTIONS PLACEHOLDER_', ligandMembraneInteractions)
		scriptTemp = scriptTemp.replace('_MOLECULAR DYNAMICS DUMP PLACEHOLDER_', 'dump			coords all custom {} {} id type x y z c_cls\ndump_modify	coords sort id'.format(self.dumpres, os.path.join(self.outdir, self.outName)))	
		scriptTemp = scriptTemp.replace('_TIMESTEP PLACEHOLDER_', 'timestep       {}'.format(self.timestep))	
		scriptTemp = scriptTemp.replace('_RUNTIME PLACEHOLDER_', 'run            {}'.format(self.run))

		with open(simData, 'w') as dataFile:
			dataFile.write(dataTemp)

		with open(simScript, 'w') as scriptFile:
			scriptFile.write(scriptTemp)
			
		print "saved files: "+str(simData)+", "+str(simScript)

	def postProcessOutput(self,outPath):
		#this function adds the ligand strengths as a column in the output file
		headerLength = 8
		lowestLigandNumber = 3
		if not os.path.exists(outPath):
			return
		with open(outPath) as f:
			lines = f.readlines()
		lines = [x.strip() for x in lines]
		skipLines = 0
		for i in range(len(lines)):
			if str('ITEM: TIMESTEP') in lines[i]:
				skipLines=headerLength
			if skipLines>0:
				skipLines-=1
				continue
			if str('ITEM: ATOMS') in lines[i]:
				lines[i] += ' affinity'
			else:
				cols = lines[i].split(' ')
				if len(cols) < 2:
					continue
				atomType = misctools.toInt(cols[1])
				if atomType<lowestLigandNumber:
					#the atom is the vehicle or a part of the membrane
					lines[i] += ' 0.0'
					#god, this is a mess...
				else:
					ligandType = atomType - lowestLigandNumber
					lines[i] += ' ' + str(self.protein.ligands[ligandType].eps)

		#add the new lines back in
		lines = [x+'\n' for x in lines]

		#write file as an xyza file
		with open(outPath+'a','w') as f:
			for l in lines:
				f.write(l)




	def deleteFiles(self):
		os.remove(os.path.join(self.filedir, self.scriptName))
		os.remove(os.path.join(self.filedir, self.dataName))
		print "removed files: "+str(os.path.join(self.filedir, self.scriptName))+", "+str(os.path.join(self.filedir, self.dataName))
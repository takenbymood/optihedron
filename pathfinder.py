import os
import argparse
from tools import analysistools as atools

import pandas as pd
from scipy.spatial import distance_matrix

parser = argparse.ArgumentParser(description='')


parser.add_argument('-i','--input', default=None, type=str, 
					help='Path of input or directory containing the input files. These must have the .xyza extension')
parser.add_argument('-o','--out', default='', type=str, 
					help='Output directory path')



def main():

	args = parser.parse_args()
	inPath = args.input
	outPath = args.out

	if not os.path.isdir(outPath):
		os.makedirs(outPath)

	jitterPath = os.path.join(outPath,'jitter/')
	matPath = os.path.join(outPath,'matrices/')
	capPath = os.path.join(outPath,'captures/')
	freePath = os.path.join(outPath,'frees/')

	if not os.path.isdir(jitterPath):
		os.makedirs(jitterPath)
	if not os.path.isdir(matPath):
		os.makedirs(matPath)
	if not os.path.isdir(capPath):
		os.makedirs(capPath)
	if not os.path.isdir(freePath):
		os.makedirs(freePath)

	if os.path.isdir(inPath):
		xyzas = [os.path.join(inPath,f) for f in os.listdir(inPath) if '.' in f and f.split('.')[-1] == 'xyza']
		
		for filePath in xyzas:
			name = filePath.split('.')[0]
			f = atools.readXYZA(filePath)

			ligands = [a[0] for a in f['atoms'] if a[1]>2]

			# get core location
			core = [a[0] for a in f['atoms'] if a[1]==2][-1]
			cxyz = [(a['x'],a['y'],a['z']) for a in f['steps'][0]['data'] if a['id'] == core][-1]

			# get relative positions of ligands
			polPos = [(a['id'],atools.crt2SphPol((a['x']-cxyz[0],a['y']-cxyz[1],a['z']-cxyz[2]))[1:] )for a in f['steps'][0]['data'] if a['id'] in ligands and a['affinity'] > 0.0]

			greatArcs = {}
			for p in polPos:
				greatArcs[p[0]] = {}
				for q in polPos:
					if p[0] == q[0]:
						greatArcs[p[0]][q[0]] = 0.0
					elif q[0] in greatArcs:
						greatArcs[p[0]][q[0]] = greatArcs[q[0]][p[0]]
					else:
						greatArcs[p[0]][q[0]] = atools.greatArcDist(p[1],q[1],rad=4.0)

			

			capSteps = []

			for step in f['steps']:
				capSteps.append((step['t'],[a['id']for a in step['data'] if a['type'] > 2 and a['c_cls']==1 and a['affinity'] > 0.0]))

			capData = []
			freeData = []
			jitterData = []

			for i in range(len(capSteps[1:])):
				captured = list(set(capSteps[i+1][-1]) - set(capSteps[i][-1]))
				freed = list(set(capSteps[i][-1]) - set(capSteps[i+1][-1]))
				if len(captured) > 0 or len(freed) > 0:
					for c in captured:
						capData.append(str(capSteps[i][0])+','+str(c))
					for r in freed:
						freeData.append(str(capSteps[i][0])+','+str(r))
				jitterData.append(str(capSteps[i][0])+','+str(len(captured))+','+str(len(freed)))
			
			distMatr = []
			for k1,v1 in greatArcs.iteritems():
				for k2, v2 in v1.iteritems():
					distMatr.append(str(k1)+','+str(k2)+','+str(v2).replace('[','').replace(']',''))

			jitterFile = os.path.join(jitterPath,name+"_jitter.csv")

			capFile = os.path.join(capPath,name+"_cap.csv")

			freeFile = os.path.join(freePath,name+"_free.csv")

			matFile = os.path.join(matPath,name+"_mat.csv")

			with open(jitterFile, 'w') as file:
				file.write("timestep,captured,freed\n")
				for item in jitterData:
					file.write("%s\n" % item)

			with open(capFile, 'w') as file:
				file.write("timestep,cap_id\n")
				for item in capData:
					file.write("%s\n" % item)

			with open(freeFile, 'w') as file:
				file.write("timestep,free_id\n")
				for item in freeData:
					file.write("%s\n" % item)

			with open(matFile, 'w') as file0210:
				file.write("i,j,great_arc_distance\n")
				for item in distMatr:
					file.write("%s\n" % item)


		

if __name__ == "__main__":
	main()
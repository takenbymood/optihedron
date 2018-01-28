import os
import shutil

baseDir = os.path.dirname(__file__)
scratch = os.path.join(baseDir,'stitchmembrane.tmp')
template = os.path.join(baseDir, 'basedata.template')
filledTemplate = os.path.join(baseDir, 'data.template')
coordsHeaderSize = 2
dipoleHeaderSize = 9

def fillTemplate(template, filledTemplate, placeHolder, filledPlaceHolder):
	with open(template) as templateFile:
		with open(filledTemplate, 'w') as filledTemplateFile:
			for line in templateFile:
				filledTemplateFile.write(line.replace(placeHolder, filledPlaceHolder))
	os.remove(template)
	os.rename(filledTemplate, template)

with open(os.path.join(baseDir, 'relaxmembrane.xyz')) as f:
	xyz = f.readlines()[coordsHeaderSize:]
with open(os.path.join(baseDir, 'relaxmembrane.mu')) as f:
	mu = f.readlines()[dipoleHeaderSize:]

membranePosition = ''
for index, (xyz_i, mu_i) in enumerate(zip(xyz, mu),1):	
	xyzpieces = xyz_i.strip().split(' ')
	mupieces = mu_i.strip().split(' ')
	membranePosition += '{0} {1} {2} {3} {4} 1 1 0'.format(index, *xyzpieces)
	membranePosition += ' {0} {1} {2}\n'.format(*mupieces)
membranePosition = membranePosition[:-2] #remove trailing newline
shutil.copyfile(template, filledTemplate)
fillTemplate(filledTemplate, scratch, '_MEMBRANE POSITIONS PLACEHOLDER_', membranePosition)
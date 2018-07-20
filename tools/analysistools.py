import numpy as np

def dropChildren(data, parentKey, childKeys, silent=True):
	if not silent:
		print('removing {} from {}!'.format(childKeys, parentKey))
	for parentItem in data[parentKey]:
		for childKey in childKeys:
			del parentItem[childKey]
	return data

def dropParentsConditional(data, parentKey, condition):
	data[parentKey] = [i for i in data[parentKey] if condition(i)]
	return data

def dropParents(data, parentKey):
	del data[parentKey]
	return data

def cleanLigands(data):
    for ind in data['individuals']:
        ind['phenome'].particle.ligands = [lig for lig in ind['phenome'].particle.ligands if lig.eps != 0]
    return data

def cleanFits(data):
    for ind in data['individuals']:
        if ind['fitness'] > 400.0:
            ind['fitness'] = 400.0 + 600.0*(1.0-((np.sum([lig.eps for lig in ind['phenome'].particle.ligands]))/(15.0*72.0)))
    return data

def sortIndLigandCount(data):
    data['individuals'].sort(key = lambda ind: len(ind['phenome'].particle.ligands))
    return data

def sortIndGens(data):
    data['individuals'].sort(key = lambda ind: ind['gen'])
    return data

def sortIndFits(data):
    data['individuals'].sort(key = lambda ind: ind['fitness'])
    return data

def load(DBs, dbConnect, clean = None, drop = None, sort = None):
	data = []

	for DB in DBs:
		dbConn = dbConnect.DatabaseConnection(DB)
		datum = dbConn.loadSession('1')
		
		if clean:
			for task in clean:
				datum = task(datum)
		if drop:
			for task in drop:
				if isinstance(task, basestring):
					del datum[task]
				else:
					datum = task(datum)
		# if sort:
		# 	for task in sort:
		# 		datum = task(datum)

		data.append(datum)
		dbConn.close()
		
	return data
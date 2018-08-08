import json
import generators

epstotal = 150
numligandtiles = 72
epsmin = 14

initPop = []
pop = 78
demes = 2

for d in range(demes):
    popTmp = []
    for p in range(pop):
        popTmp.append(generators.fixedAffinityGenomes(epstotal, epsmin, numligandtiles))        
    initPop.append(popTmp)


initParams = {'seed': 123456, 'partialpacking': False, 'pop': pop, 'demes': demes, 'epsmin': epsmin, 'init_pop': initPop}


initFileName = '{}-{}.json'.format(epstotal,initParams['epsmin'])
with open(initFileName, 'w') as initFile:
    json.dump(initParams, initFile)
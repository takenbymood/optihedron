import json
import numpy as np

epstotal = 150
numligandtiles = 72
epsmin = 14

def guessGenome(epstotal,epsmin):    
    numLigandsExact = float(epstotal)/float(epsmin)
    numLigandsApprox = int(numLigandsExact)
    numLigandsRemainder = numLigandsExact - numLigandsApprox
    numLigandsInit = numLigandsApprox + np.random.choice([0,1],p=[1.0-numLigandsRemainder,numLigandsRemainder])

    cleanGenome = [0] * numligandtiles    
    replacementPos = np.random.choice(len(cleanGenome), size=numLigandsInit, replace=False)
    for i in replacementPos:
        cleanGenome[i] = 1

    return cleanGenome

initPop = []
pop = 78
demes = 2

for d in range(demes):
    popTmp = []
    for p in range(pop):
        popTmp.append(guessGenome(epstotal, epsmin))
    initPop.append(popTmp)


initParams = {'seed': 123456, 'partialpacking': False, 'pop': pop, 'demes': demes, 'epsmin': epsmin, 'init_pop': initPop}


initFileName = '{}-{}.json'.format(epstotal,initParams['epsmin'])
with open(initFileName, 'w') as initFile:
    json.dump(initParams, initFile)
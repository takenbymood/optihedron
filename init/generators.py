import numpy as np

def fixedAffinityGenomes(epstotal,epsmin,numligandtiles):    
    numLigandsExact = float(epstotal)/float(epsmin)
    numLigandsApprox = int(numLigandsExact)
    numLigandsRemainder = numLigandsExact - numLigandsApprox
    numLigandsInit = numLigandsApprox + np.random.choice([0,1],p=[1.0-numLigandsRemainder,numLigandsRemainder])

    cleanGenome = [0] * numligandtiles    
    replacementPos = np.random.choice(len(cleanGenome), size=numLigandsInit, replace=False)
    for i in replacementPos:
        cleanGenome[i] = 1

    return cleanGenome

def fixedLigandsGenomes(numligs,epsmin,numligandtiles):    
    numLigandsInit = numligs

    cleanGenome = [0] * numligandtiles    
    replacementPos = np.random.choice(len(cleanGenome), size=numLigandsInit, replace=False)
    for i in replacementPos:
        cleanGenome[i] = 1

    return cleanGenome
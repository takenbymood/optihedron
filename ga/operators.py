import random
import numpy as np
from deap import tools

def defaultMut(individual, MINPDB):    
    return tools.mutFlipBit(individual, MINPDB)

def fixedActivationMut(individual, MINPDB):   
    return tools.mutShuffleIndexes(individual, MINPDB) 

def defaultGeneWiseTwoPoint(ind1, ind2, GENES, GENESIZE):              
    cxpoint1 = random.randint(1, GENES)*GENESIZE
    cxpoint2 = random.randint(1, GENES-1)*GENESIZE
    if cxpoint2 >= cxpoint1:
        cxpoint2 += GENESIZE
    else: # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
        
    return ind1, ind2

def fixActivation(ind, targetActivations):

    misActivatedCount = np.sum(ind) - targetActivations    
    misActivatedState = 1 if misActivatedCount > 0 else 0    
    targetActivationState = 1 if misActivatedState == 0 else 0    
    misActivatedPos = [idx for idx, val in enumerate(ind) if val == misActivatedState]    
    fixActivationPos = np.random.choice(misActivatedPos, size=abs(misActivatedCount), replace=False)    

    for pos in fixActivationPos:
        ind[pos] = targetActivationState

    return ind

def fixedActivationGeneWiseTwoPoint(ind1, ind2, GENES, GENESIZE):    
    activates1 = np.sum(ind1)
    activates2 = np.sum(ind2)

    ind1, ind2 = defaultGeneWiseTwoPoint(ind1, ind2, GENES, GENESIZE)
    ind1, ind2 = fixActivation(ind1, activates1), fixActivation(ind2, activates2)    

    return ind1, ind2
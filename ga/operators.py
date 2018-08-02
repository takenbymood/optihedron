import random
from deap import tools

def defaultMut(individual, MINPDB):    
    return tools.mutFlipBit(individual, MINPDB)

def defaultGeneWiseTwoPoint(ind1,ind2, GENES, GENESIZE):            
    cxpoint1 = random.randint(1, GENES)*GENESIZE
    cxpoint2 = random.randint(1, GENES-1)*GENESIZE
    if cxpoint2 >= cxpoint1:
        cxpoint2 += GENESIZE
    else: # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
        
    return ind1, ind2
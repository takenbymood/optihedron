import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import copy
import time
import pickle

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
		if sort:
			for task in sort:
				datum = task(datum)

		data.append(datum)
		dbConn.close()
		
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

def sameLigands(ligA, ligB):
	if ligA.__dict__ == ligB.__dict__:
		return True
	else:
		return False

def greatArcDist(Ang1, Ang2, rad=4):
    #Ang = (PolarAng,AziAng)
    #https://math.stackexchange.com/questions/231221/great-arc-distance-between-two-points-on-a-unit-sphere
    arcDist=rad*(np.arccos((np.cos(Ang1[0])*np.cos(Ang2[0]))+((np.sin(Ang1[0])*np.sin(Ang2[0]))*(np.cos(Ang1[1]-Ang2[1])))))
    return arcDist

def tagLigands(data):
	data['individuals'] = [classifyLigands(ind) for ind in data['individuals']]
	return data		

def classifyLigands(ind):
    cutoff = 2.0

    if 'ligtags' in ind:
		return ind['ligtags']
    else:
        ligands = ind['phenome'].particle.ligands
        spotty = 0
        liney = 0
        patchy = 0
        identity = {}

        for ligCursor in ligands:	
            identity[ligCursor] = {}
            NNcount = 0
            NNlist = []

            for ligOther in ligands:
                if not sameLigands(ligCursor, ligOther):
                    arcDist = greatArcDist((ligCursor.polAng,ligCursor.aziAng),(ligOther.polAng,ligOther.aziAng)) 
                    if arcDist <= cutoff:
                        NNcount += 1
                        NNlist.append(ligOther)

            identity[ligCursor]['NNcount'] = NNcount
            identity[ligCursor]['NNlist'] = NNlist

            if NNcount == 3: #if 3 NN
                arcDist12 = greatArcDist((NNlist[0].polAng,NNlist[0].aziAng),(NNlist[1].polAng,NNlist[1].aziAng))
                arcDist13 = greatArcDist((NNlist[0].polAng,NNlist[0].aziAng),(NNlist[2].polAng,NNlist[2].aziAng))
                arcDist23 = greatArcDist((NNlist[1].polAng,NNlist[1].aziAng),(NNlist[2].polAng,NNlist[2].aziAng))

                if not arcDist12 <= 2.0 and not arcDist13 <= 2.0 and not arcDist23 <= 2.0: #and noone of them are ea other's NN                    
                    liney += 1
                    identity[ligCursor]['character'] = 'liney'

            if NNcount == 2:#if 2 NN                
                arcDist = greatArcDist((NNlist[0].polAng,NNlist[0].aziAng),(NNlist[1].polAng,NNlist[1].aziAng))
                if not arcDist <= 2.0:#and the NNs are not each other's NN
                    liney += 1
                    identity[ligCursor]['character'] = 'liney'
          
            if NNcount == 1:#if 1 NN
                liney += 1
                identity[ligCursor]['character'] = 'liney'

            if NNcount == 0: #if 0 NN
                spotty += 1
                identity[ligCursor]['character'] = 'spotty'
            
            #all other ligands contribute to overall patchy-ness            
            if 'character' not in identity[ligCursor]:
                identity[ligCursor]['character'] = 'patchy'

        patchy = len(ligands) - spotty - liney
        ind['ligtags'] = patchy, spotty, liney, identity

        return patchy, spotty, liney, identity

def clusterLineyLigands(ligands, silent=True):
    #print(greatArcDist((0,0), ((3.14/15.0),0))) #bad at small distances, correct should be 0.8377 but returns 0.8373, willing to accept +- 0.001
    #ligands = [(ligand,ligandTags)]
    if not silent:
        startTime = time.time()
    
    ligandsTmp = copy.deepcopy(ligands)

    clusters = []
    while ligandsTmp:            
        nextSeedQueue = [ligandsTmp[0]]    

        clusterTmp = []
        clusterTmp.append(ligandsTmp[0])
        del ligandsTmp[0]
        while nextSeedQueue:
            seed = nextSeedQueue.pop()
            rLigands = []
            for ligand in ligandsTmp:
                NNlist = ligand[1]['NNlist']
                existingCopies = 0
                for NN in NNlist:
                    if sameLigands(seed[0], NN):
                        existingCopies += 1
                if existingCopies == 1:                                    
                    nextSeedQueue.append(ligand)
                    clusterTmp.append(ligand)
                    rLigands.append(ligands)
                elif existingCopies > 1:
                    raise ValueError
            for rLig in rLigands:
                ligandsTmp.remove(rLig)
        clusters.append(clusterTmp)        

    if not silent:
        print('ligands: {}'.format(len(ligands)))
        print('clusters: {}'.format(len(clusters)))
        print('clustered in {}s'.format(time.time()-startTime))

    return clusters

def groupLineyLigands(ind):
    if 'lineyTags' in ind:
        return ind['lineyTags']
    else:
        _, _, _, identity = classifyLigands(ind)
        lineyLigands = []

        for ligand, ligandTags in identity.items():
            if ligandTags['character'] == 'liney':
                lineyLigands.append((ligand,ligandTags)) 
        lineyClusters = clusterLineyLigands(lineyLigands)
        ind['lineyTags'] = lineyClusters
        return lineyClusters

def scanGen(scanData, interest, indexOffset, aggregateMode, silent=True):
    scanPlotIndex = []
    scanPlotData = []

    aggregateMethod = None
    if aggregateMode != 'MIN' and aggregateMode != 'MAX' and aggregateMode != 'AVG':
        raise ValueError    

    if not silent:
        startTime = time.time()
    for scanDatum in scanData:

        aggregateInterest = []
        aggregateInds = 0
        cursorGen = 0
        genAveragedInterest = []        

        expectedTicks = len(scanDatum['individuals'])
        actualTicks = 0

        for ind in scanDatum['individuals']:            
            if ind['gen'] != cursorGen:    
                if not silent:
                    print('analyzed gen{}'.format(cursorGen))

                if aggregateMode == 'MIN':
                    genAveragedInterest.append(((cursorGen),np.min(aggregateInterest)))
                elif aggregateMode == 'MAX':
                    genAveragedInterest.append(((cursorGen),np.max(aggregateInterest)))
                elif aggregateMode == 'AVG':
                    genAveragedInterest.append(((cursorGen),(np.sum(aggregateInterest)/float(aggregateInds))))
                
                actualTicks += aggregateInds
                cursorGen = ind['gen']                    
                aggregateInterest = [interest(ind)]        
                aggregateInds = 1
            else:                
                aggregateInterest.append(interest(ind)) 
                aggregateInds += 1

        if aggregateMode == 'MIN':
            genAveragedInterest.append(((cursorGen),np.min(aggregateInterest)))
        elif aggregateMode == 'MAX':
            genAveragedInterest.append(((cursorGen),np.max(aggregateInterest)))
        elif aggregateMode == 'AVG':
            genAveragedInterest.append(((cursorGen),(np.sum(aggregateInterest)/float(aggregateInds))))
        actualTicks += aggregateInds

        if expectedTicks != actualTicks:
            raise ValueError

        scanPlotIndex.append(indexOffset)
        scanPlotData.append(genAveragedInterest)            
        if not silent:
            print('analyzed gen{}'.format(cursorGen))
            print('analyzed index {} of scanData'.format(indexOffset))
            print('analysis took {}s'.format(time.time() - startTime))            
        indexOffset += 1
    return scanPlotIndex, scanPlotData

def scanCustom(scanData, interestKey, interestKeyLabel, tickerRange, tickerBlockSize, tickerBlockOffset, interest, indexOffset, aggregateMode, tickerInterval = 0.0, silent=True):
    scanPlotIndex = []
    scanPlotData = []    

    aggregateMethod = None
    if aggregateMode != 'MIN' and aggregateMode != 'MAX' and aggregateMode != 'AVG':
        raise ValueError    

    if not silent:
        startTime = time.time()
    for scanDatum in scanData:

        ticks = []
        for tickerBlock in range(tickerRange):
            ticks.append(tickerBlock*tickerBlockSize+tickerBlockOffset)        

        aggregateInterest = []
        aggregateInds = 0
        cursorTick = ticks.pop(0)
        tickAveragedInterest = []        

        expectedTicks = len(scanDatum['individuals'])
        actualTicks = 0        

        for ind in scanDatum['individuals']:                                   
            if interestKey(ind) - cursorTick > tickerInterval:                    
                if not silent:
                    print('analyzed tick {} of {}'.format(cursorTick, interestKeyLabel))

                if aggregateInds != 0:
                    if aggregateMode == 'MIN':
                        tickAveragedInterest.append(((cursorTick),np.min(aggregateInterest)))
                    elif aggregateMode == 'MAX':
                        tickAveragedInterest.append(((cursorTick),np.max(aggregateInterest)))
                    elif aggregateMode == 'AVG':
                        tickAveragedInterest.append(((cursorTick),(np.sum(aggregateInterest)/float(aggregateInds))))
                    
                    actualTicks += aggregateInds
                    cursorTick = ticks.pop(0)

                while interestKey(ind) - cursorTick > tickerInterval:                        
                    tickAveragedInterest.append(((cursorTick),np.nan))
                    cursorTick = ticks.pop(0)
                
                aggregateInterest = [interest(ind)]        
                aggregateInds = 1
            else:                
                aggregateInterest.append(interest(ind)) 
                aggregateInds += 1

        if aggregateMode == 'MIN':
            tickAveragedInterest.append(((cursorTick),np.min(aggregateInterest)))
        elif aggregateMode == 'MAX':
            tickAveragedInterest.append(((cursorTick),np.max(aggregateInterest)))
        elif aggregateMode == 'AVG':
            tickAveragedInterest.append(((cursorTick),(np.sum(aggregateInterest)/float(aggregateInds))))
        actualTicks += aggregateInds

        if expectedTicks != actualTicks:            
            raise ValueError

        while ticks:
            tickAveragedInterest.append(((ticks.pop(0)),np.nan))

        scanPlotIndex.append(indexOffset)
        scanPlotData.append(tickAveragedInterest)            
        if not silent:
            print('analyzed tick {} of {}'.format(cursorTick, interestKeyLabel))
            print('analyzed index {} of scanData'.format(indexOffset))
            print('analysis took {}s'.format(time.time() - startTime))            
        indexOffset += 1
    return scanPlotIndex, scanPlotData

def characterScore(ind, patchyScore = -1, spottyScore = -1, lineyScore = 1):
    patchy, spotty, liney, _ = classifyLigands(ind)
    charScore = float((patchy*patchyScore) + (spotty*spottyScore) + (liney*lineyScore))/float(patchy+spotty+liney)

    return charScore

def lineyABS(ind):
    _, _, liney, _ = classifyLigands(ind)

    return float(liney)

def spottyABS(ind):
    _, spotty, _, _ = classifyLigands(ind)

    return float(spotty)

def patchyABS(ind):
    patchy, _, _, _ = classifyLigands(ind)

    return float(patchy)

def lineyREL(ind):
    patchy, spotty, liney, _ = classifyLigands(ind)

    return float(liney)/float(patchy+spotty+liney)

def spottyREL(ind):
    patchy, spotty, liney, _ = classifyLigands(ind)

    return float(spotty)/float(patchy+spotty+liney)

def patchyREL(ind):
    patchy, spotty, liney, _ = classifyLigands(ind)

    return float(patchy)/float(patchy+spotty+liney)

def NNCountAVG(ind):
    _, _, _, identity = classifyLigands(ind)
    totalNNs = 0
    totalLigands = 0
    for ligand, ligandTags in identity.items():
        totalNNs += ligandTags['NNcount']
        totalLigands += 1
    return float(totalNNs)/float(totalLigands)

def lineyChainTags(ind, minChainLength):
    lineyClusters = groupLineyLigands(ind)    
    lineyChains = []

    for cluster in lineyClusters:
        if len(cluster) >= minChainLength:            
            lineyChains.append(cluster)

    return lineyChains    

def lineyChainCount(ind, minChainLength=3):
    lineyChains = lineyChainTags(ind, minChainLength)

    return float(len(lineyChains))

def lineyChainSizeAVG(ind, minChainLength=3):
    lineyChains = lineyChainTags(ind, minChainLength)

    if lineyChains:
        return float(np.average([len(chain) for chain in lineyChains]))
    else:
        return 0.0

def lineyChainSizeMAX(ind, minChainLength=3):
    lineyChains = lineyChainTags(ind, minChainLength)

    if lineyChains:
        return float(np.max([len(chain) for chain in lineyChains]))
    else:
        return 0.0

def lineyChainSizeMIN(ind, minChainLength=3):
    lineyChains = lineyChainTags(ind, minChainLength)

    if lineyChains:
        return float(np.min([len(chain) for chain in lineyChains]))
    else:
        return 0.0

def plotScanGen(scanData, scanLabel, scanIndices, interest, indexOffset, aggregateMode, plotName, cmap, fmt='.2g', vmin = None, vmax = None, annotate=False, silent=True, visual=True, dump=False, backup=False, dumpdir='plots'):    
    if not os.path.exists(dumpdir):
        os.mkdir(dumpdir)

    scanPlotIndex, scanPlotData = scanGen(scanData, interest, indexOffset, aggregateMode, silent=silent)

    if backup:
        pickle.dump(scanPlotIndex, open("{}-gen-scanPlotIndex.pickle".format(plotName), "wb"))
        pickle.dump(scanPlotData, open("{}-gen-scanPlotData.pickle".format(plotName), "wb"))

    plotData = np.zeros((len([i[1] for i in scanPlotData[0]]),len(scanPlotIndex)))
    annotData = np.zeros((len([i[1] for i in scanPlotData[0]]),len(scanPlotIndex)))

    cursorA = 0
    for i in range(len(scanPlotIndex)):
        scanPlotDatum = scanPlotData[cursorA]
        cursorB = 0
        for _ in [i[1] for i in scanPlotDatum]:        
            plotData[cursorB][cursorA] = scanPlotDatum[cursorB][1]
            annotData[cursorB][cursorA] = scanPlotDatum[cursorB][1]
            cursorB += 1
        cursorA += 1

    annot = annotData if annotate else False        

    ax = sns.heatmap(plotData, linewidth=1, annot=annot, cmap=cmap, vmin=vmin, vmax=vmax, fmt=fmt)
    plt.title('{}'.format(plotName))
    ax.set_xticks([i+0.5-indexOffset for i in scanIndices])
    ax.invert_yaxis()
    ax.set_xticklabels([i for i in scanIndices])
    plt.xlabel('{}'.format(scanLabel))
    plt.ylabel('generation')
    if dump:
        plt.savefig('{}/{}.png'.format(dumpdir,plotName));
    if visual:
        plt.show();

    
def plotScanCustom(scanData, scanLabel, scanIndices, interest, indexOffset, aggregateMode, plotName, cmap, interestKey, interestKeyLabel, tickerRange, tickerBlockSize, tickerBlockOffset, tickerInterval=0.0, fmt='.2g', vmin = None, vmax = None, annotate=False, linecolor='black', silent=True, visual=True, dump=False, backup=False, dumpdir='plots'):    
    if not os.path.exists(dumpdir):
        os.mkdir(dumpdir)

    scanPlotIndex, scanPlotData = scanCustom(scanData, interestKey, interestKeyLabel, tickerRange, tickerBlockSize, tickerBlockOffset, interest, indexOffset, aggregateMode, tickerInterval=tickerInterval,silent=silent)

    if backup:
        pickle.dump(scanPlotIndex, open("{}-{}-scanPlotIndex.pickle".format(plotName,interestKeyLabel), "wb"))
        pickle.dump(scanPlotData, open("{}-{}-scanPlotData.pickle".format(plotName,interestKeyLabel), "wb"))

    plotData = np.zeros((len([i[1] for i in scanPlotData[0]]),len(scanPlotIndex)))
    annotData = np.zeros((len([i[1] for i in scanPlotData[0]]),len(scanPlotIndex)))

    cursorA = 0
    for i in range(len(scanPlotIndex)):
        scanPlotDatum = scanPlotData[cursorA]
        cursorB = 0
        for _ in [i[1] for i in scanPlotDatum]:        
            plotData[cursorB][cursorA] = scanPlotDatum[cursorB][1]
            annotData[cursorB][cursorA] = scanPlotDatum[cursorB][1]
            cursorB += 1
        cursorA += 1

    annot = annotData if annotate else False        

    ax = sns.heatmap(plotData, linewidth=1, annot=annot, cmap=cmap, vmin=vmin, vmax=vmax, fmt=fmt, linecolor=linecolor)
    plt.title('{}'.format(plotName))
    ax.set_xticks([i+0.5-indexOffset for i in scanIndices])
    ax.set_xticklabels([i for i in scanIndices])
    ax.invert_yaxis()
    yticklabels = [0]
    for tickerBlock in range(tickerRange):
        yticklabels.append(tickerBlock*tickerBlockSize+tickerBlockOffset)
    ax.set_yticks([i for i in range(0,(tickerRange+1))])
    ax.set_yticklabels(yticklabels)
    ax.set_facecolor('#F5F5F5')
    plt.xlabel('{}'.format(scanLabel))
    plt.ylabel('{}'.format(interestKeyLabel))
    if dump:
        plt.savefig('{}/{}.png'.format(dumpdir,plotName));
    if visual:
        plt.show();

def measureLigandContact(xyzaFile, headersize=9, xyzallsize=2973, timestepinterval=100, intrange=1.8, silent=True):
    with open(xyzaFile) as f:
        content = f.readlines()
    content = [x.strip() for x in content] 

    steps = []
    timestep = 0
    while content:
        #remove header    
        xyz_iter = content[headersize:headersize+xyzallsize]
        steps.append((timestep, xyz_iter))
        content = content[headersize+xyzallsize:]
        timestep += timestepinterval

    contactData = []
    progress = 1
    for step in steps:
        if not silent:
            print('{}/{}'.format(progress, len(steps)))
        progress+=1
        time, coords = step
        coordsCLEANED = []
        for coord in coords:
            coordsCLEANED.append(coord.split())
        ligands = [i for i in coordsCLEANED if float(i[-1]) > 0.0]
        mems = [i for i in coordsCLEANED if float(i[1]) == 1.0]
        
        cLIG = []
        cMEM = []
        #convert strings to ints...    
        for ligand in ligands:
            ligand = [float(i) for i in ligand]
            cLIG.append(ligand)
        for mem in mems:
            mem = [float(i) for i in mem]
            cMEM.append(mem)
        ligands = cLIG
        mems = cMEM
        
        #count number of ligands that are within range
        ligandContact = 0
        for ligand in ligands:
            contact = False
            for mem in mems:
                if not contact:
                    lx, ly, lz = ligand[2], ligand[3], ligand[4]
                    mx, my, mz = mem[2], mem[3], mem[4]
                    dist = np.sqrt( (lx - mx)**2 + (ly - my)**2 + (lz - mz)**2 )
                    if dist <= intrange:
                        contact = True
            if contact:
                ligandContact += 1
        if ligandContact > len(ligands):
            raise ValueError
        contactData.append((time, ligandContact))

    return contactData, len(ligands)
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import copy
import time
import pickle
import networkx
import holoviews as hv



def buildNanoParticleFromNetwork(G,weight,radius=4,sig=1):
    particle = nanoparticle.NanoParticle()
    for n,w in G.nodes(data=True):
        particle.ligands.append(nanoparticle.Ligand(weight,sig,radius,w['pol'],w['azi']))
    return particle


def buildLigandNetwork(ligands, silent=True):
    if not silent:
        startTime = time.time()
    G=networkx.Graph()
    
    nIndex = 1
    for i in ligands:
        if i.eps > 0.0:
            G.add_node(nIndex,weight=i.eps,polAng=i.polAng,aziAng=i.aziAng)
        nIndex += 1
    
    iIndex = 1
    for i in ligands:
        jIndex = 1
        for j in ligands:
            if i < j:
                cartDist = 1.0/greatArcDist((i.polAng,i.aziAng),(j.polAng,j.aziAng))
                #affDist = abs(i.eps - j.eps)
                if i.eps > 0.0 and j.eps > 0.0:
                    G.add_edge(iIndex, jIndex, weight=cartDist)
            jIndex += 1
        iIndex += 1
    return G

def pruneNetwork(G,pruning):
    prunes = []
    GP = copy.deepcopy(G)
    maxW = 0
    for e in GP.edges:
        w = GP.get_edge_data(*e)['weight']
        if(w<=pruning):
            prunes.append(e)

    GP.remove_edges_from(prunes)
    return GP


def buildPrunedLigandNetwork(ind,pruning):
    G = buildLigandNetwork(ind.phenomePickle.particle.ligands)
    return pruneNetwork(G,pruning)

def buildPrunedNetworkView(ind,pruning):
    G = buildPrunedLigandNetwork(ind,pruning)
    padding = dict(x=(-1.2, 1.2), y=(-1.2, 1.2))
    return hv.Graph.from_networkx(G, networkx.layout.spring_layout).opts(plot=dict(color_index='weight')).redim.range(**padding)

def buildNetworkView(G):
    padding = dict(x=(-1.2, 1.2), y=(-1.2, 1.2))
    return hv.Graph.from_networkx(G, networkx.layout.spring_layout).opts(plot=dict(color_index='weight')).redim.range(**padding)

def hammingDistance(s1,s2):
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def formFactor(ligands):
    
    spottyLigands = 0    
    lineyLigands = 0    
    patchyLigands = 0
    
    minDist = 2.0
    smallestDist = 100

    #counting liney ligands
    nZero = 0
    for j in ligands:
        if j.eps == 0.0:
            nZero+=1
            continue
        NNlist = []
        for k in ligands:        
            if j != k and k.eps > 0.0:
                dist = greatArcDist((j.polAng,j.aziAng),(k.polAng,k.aziAng))
                if dist <= minDist:
                    NNlist.append(k)
                if dist < smallestDist:
                    smallestDist = dist
        NNcount = len(NNlist)
        if NNcount == 2:#if it has 2 NNs 
            dist = greatArcDist((NNlist[0].polAng,NNlist[0].aziAng),(NNlist[1].polAng,NNlist[1].aziAng))
            if not dist <= minDist:#and the NNs are not each other's NN
                lineyLigands += 1                              
        elif NNcount == 1:
            lineyLigands += 1 
        elif NNcount == 0:
            spottyLigands += 1
    
    #counting patchy ligands
    totalLigands = float(len(ligands)) - nZero
    patchyLigands = totalLigands - spottyLigands - lineyLigands
    
    return [float(patchyLigands)/totalLigands, float(lineyLigands)/totalLigands, float(spottyLigands)/totalLigands]

def buildNetworkList(individuals):
    networks = []
    for i in individuals:
        n = buildLigandNetwork(i.phenomePickle.particle.ligands)
        networks.append((n,i))
    return networks

def pruneNetworkList(networks,pruning):
    prunedNetworks = []
    for n in networks:
        prunedNetworks.append((pruneNetwork(n[0],pruning),n[1]))
    return prunedNetworks

def buildPrunedNetworkList(individuals,pruning):
    return pruneNetworkList(buildNetworkList(individuals),pruning)

def gnmw(G):
    return(len(G.nodes()),len(G.edges()),np.sum([w['weight'] for (u,v,w) in G.edges(data=True)]))

def gnmRandomWeightedGraph(nodes,edges,weight):
    G = networkx.gnm_random_graph(nodes,edges)
    totalWeight = 0
    for (u,v,w) in G.edges(data=True):
        w['weight'] = random.uniform(0,10)
        totalWeight += w['weight']
    
    normFactor = weight/totalWeight
    
    for (u,v,w) in G.edges(data=True):
        w['weight'] = w['weight']*normFactor
        
    return G

def transparent_cmap(cmap, N=255):
    "Copy colormap and set alpha values"

    mycmap = cmap
    mycmap._init()
    mycmap._lut[:,-1] = np.linspace(0, 0.8, N+4)
    return mycmap

def smallWorldNess(G):
    gnmwG = gnmw(G)
    RG = gnmRandomWeightedGraph(gnmwG[0],gnmwG[1],gnmwG[2])
    pLengthG = networkx.average_shortest_path_length(G,weight='weight')
    pLengthRG = networkx.average_shortest_path_length(RG,weight='weight')
    clusteringG = networkx.average_clustering(G,weight='weight')
    clusteringRG = networkx.average_clustering(RG,weight='weight')
    SW = (clusteringG/clusteringRG)/(pLengthG/pLengthRG)
    return SW


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
                    rLigands.append(ligand)
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
    if 'lineyclustertags' in ind:
        return ind['lineyclustertags']
    else:
        _, _, _, identity = classifyLigands(ind)
        lineyLigands = []

        for ligand, ligandTags in identity.items():
            if ligandTags['character'] == 'liney':
                lineyLigands.append((ligand,ligandTags)) 
        lineyClusters = clusterLineyLigands(lineyLigands)
        ind['lineyclustertags'] = lineyClusters
        return lineyClusters

def scanGen(scanData, interest, indexOffset, aggregateMode, silent=True):
    scanPlotIndex = []
    scanPlotData = []

    aggregateMethod = None
    if aggregateMode != 'MIN' and aggregateMode != 'MAX' and aggregateMode != 'AVG' and aggregateMode != 'POP':
        raise ValueError    

    if not silent:
        startTime = time.time()
    for scanDatum in scanData:

        aggregateInterest = []
        aggregateBucket = []
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
                elif aggregateMode == 'POP':
                    genAveragedInterest.append(((cursorGen),interest(aggregateBucket)))
                
                actualTicks += aggregateInds
                cursorGen = ind['gen']                    
                if aggregateMode == 'MIN' or aggregateMode == 'MAX' or aggregateMode == 'AVG':
                    aggregateInterest = [interest(ind)]
                elif aggregateMode == 'POP':
                    aggregateBucket = [ind]
                aggregateInds = 1
            else:                
                if aggregateMode == 'MIN' or aggregateMode == 'MAX' or aggregateMode == 'AVG':
                    aggregateInterest.append(interest(ind)) 
                elif aggregateMode == 'POP':
                    aggregateBucket.append(ind)
                aggregateInds += 1

        if aggregateMode == 'MIN':
            genAveragedInterest.append(((cursorGen),np.min(aggregateInterest)))
        elif aggregateMode == 'MAX':
            genAveragedInterest.append(((cursorGen),np.max(aggregateInterest)))
        elif aggregateMode == 'AVG':
            genAveragedInterest.append(((cursorGen),(np.sum(aggregateInterest)/float(aggregateInds))))
        elif aggregateMode == 'POP':
            genAveragedInterest.append(((cursorGen),interest(aggregateBucket)))
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
    if aggregateMode != 'MIN' and aggregateMode != 'MAX' and aggregateMode != 'AVG' and aggregateMode != 'POP':
        raise ValueError    

    if not silent:
        startTime = time.time()
    for scanDatum in scanData:

        ticks = []
        for tickerBlock in range(tickerRange):
            ticks.append(tickerBlock*tickerBlockSize+tickerBlockOffset)        

        aggregateInterest = []
        aggregateBucket = []
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
                    elif aggregateMode == 'POP':
                        tickAveragedInterest.append(((cursorTick),interest(aggregateBucket)))
                    
                    actualTicks += aggregateInds
                    cursorTick = ticks.pop(0)

                while interestKey(ind) - cursorTick > tickerInterval:                        
                    tickAveragedInterest.append(((cursorTick),np.nan))
                    cursorTick = ticks.pop(0)
                
                if aggregateMode == 'MIN' or aggregateMode == 'MAX' or aggregateMode == 'AVG':
                    aggregateInterest = [interest(ind)]
                elif aggregateMode == 'POP':
                    aggregateBucket = [ind]
                aggregateInds = 1
            else:
                if aggregateMode == 'MIN' or aggregateMode == 'MAX' or aggregateMode == 'AVG':
                    aggregateInterest.append(interest(ind)) 
                elif aggregateMode == 'POP':
                    aggregateBucket.append(ind)                                
                aggregateInds += 1

        if aggregateMode == 'MIN':
            tickAveragedInterest.append(((cursorTick),np.min(aggregateInterest)))
        elif aggregateMode == 'MAX':
            tickAveragedInterest.append(((cursorTick),np.max(aggregateInterest)))
        elif aggregateMode == 'AVG':
            tickAveragedInterest.append(((cursorTick),(np.sum(aggregateInterest)/float(aggregateInds))))
        elif aggregateMode == 'POP':
            tickAveragedInterest.append(((cursorTick),interest(aggregateBucket)))
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

def lineyLineTags(ind, minLineLength):      
    if 'lineylinetags' in ind:
        return ind['lineylinetags']
    else:        
        lineyClusters = groupLineyLigands(ind)
        lineyLigands = []
        for cluster in lineyClusters:
            for ligand, ligandTags in cluster:
                lineyLigands.append(ligand)
        _, _, _, identity = classifyLigands(ind)
        lineyLines = []
                
        for cluster in lineyClusters:            
            if len(cluster) >= minLineLength:
                G = networkx.Graph()
                for lineyLig in cluster:
                    G.add_node(lineyLig[0])                    
                for lineyLig in cluster:
                    for NN in lineyLig[1]['NNlist']:
                        if NN in lineyLigands:                        
                            G.add_edge(lineyLig[0], NN)                        

                dists = []
                for nodeA in G.nodes():
                    tmpdists = []
                    for nodeB in G.nodes():
                        if not sameLigands(nodeA, nodeB):
                            for pathl in [len(path) for path in networkx.all_simple_paths(G, nodeA, nodeB)]:
                                tmpdists.append(pathl)                        
                    dists.append(np.max(tmpdists))                                            
                lineyLines.append(np.max(dists))

        lineyLines = [line for line in lineyLines if line >= minLineLength]

        ind['lineylinetags'] = lineyLines

        return lineyLines

def lineyLineCount(ind, minLineLength=3):
    lineyLines = lineyLineTags(ind, minLineLength)

    return float(len(lineyLines))

def lineyLineSizeAVG(ind, minLineLength=3):
    lineyLines = lineyLineTags(ind, minLineLength)

    if lineyLines:
        return float(np.average(lineyLines))
    else:
        return 0.0

def lineyLineSizeMAX(ind, minLineLength=3):
    lineyLines = lineyLineTags(ind, minLineLength)

    if lineyLines:
        return float(np.max(lineyLines))
    else:
        return 0.0

def lineyLineSizeMIN(ind, minLineLength=3):
    lineyLines = lineyLineTags(ind, minLineLength)

    if lineyLines:
        return float(np.min(lineyLines))
    else:
        return 0.0

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

def geneticDiversity(inds):
    genomes = {}    
    for ind in inds:
        genome = tuple(ind['genome'])
        if genome not in genomes:
            genomes[genome] = 1
        else:
            genomes[genome] += 1

    if len(inds) != np.sum(genomes.values()):        
        raise ValueError
    return genomes

def genomeConvergence(inds):
    genomes = geneticDiversity(inds)
    return 1.0 - float(len(genomes.keys()))/float(len(inds))  

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

def smoothLigandContact(timeWindow, ligandContactWindow, windowSize, mode):
    window = np.ones(windowSize)/windowSize
    ligandContactWindowSmooth = np.convolve(ligandContactWindow, window, mode=mode)
    timeWindow = timeWindow[windowSize:][:-windowSize]
    ligandContactWindowSmooth = ligandContactWindowSmooth[windowSize:][:-windowSize]
    return timeWindow, ligandContactWindowSmooth

def jitterLigandContact(contactData, windowSize, mode):
    time = [i[0] for i in contactData]
    ligandContact = [i[1] for i in contactData]

    timeJitter = time[windowSize:][:-windowSize]
    ligandContactInterest = ligandContact[windowSize:][:-windowSize]
    _, ligandContactSmooth = smoothLigandContact(time, ligandContact, windowSize, mode)
    
    jitter = []
    for ligandContactInterest_i, ligandContactSmooth_i in zip(ligandContactInterest, ligandContactSmooth):
        jitter.append(ligandContactInterest_i - ligandContactSmooth_i)

    return timeJitter, jitter

def jitterLigandContactSUM(contactData, windowSize=10, mode='same', jitTolerance=0.1):
    _, jitter = jitterLigandContact(contactData, windowSize, mode)
    return np.sum([abs(jit) for jit in jitter if abs(jit) >= jitTolerance])

def jitterLigandContactMAX(contactData, windowSize=10, mode='same', jitTolerance=0.1):
    _, jitter = jitterLigandContact(contactData, windowSize, mode)
    return np.max([abs(jit) for jit in jitter if abs(jit) >= jitTolerance])

def digitalJitter(jitter, jitTolerance):
    jitDigital = []
    for jit in jitter:
        if abs(jit) >= jitTolerance:
            jitDigital.append(1)
        else:
            jitDigital.append(0)
    return jitDigital

def jitterLigandContactLIFETIME(contactData, windowSize=10, mode='same', jitTolerance=0.1):
    _, jitter = jitterLigandContact(contactData, windowSize, mode)
    jitter = digitalJitter(jitter, jitTolerance)

    totalLifeTime = 0    

    jitPacketALL = []
    jitPacket = []
    for jit in jitter:
        if jit == 1:
            jitPacket.append(1)
            totalLifeTime += 1
        elif jit == 0:
            if jitPacket:
                jitPacketALL.append(jitPacket)
                jitPacket = []
        else:
            raise ValueError
    jitPacketALL.append(jitPacket)

    averageLifeTime = np.average([len(jitPack) for jitPack in jitPacketALL])
    maxLifeTime = np.max([len(jitPack) for jitPack in jitPacketALL])
    minLifeTime = np.min([len(jitPack) for jitPack in jitPacketALL])

    return totalLifeTime, averageLifeTime, maxLifeTime, minLifeTime

def metastableLIFETIME(contactData):    
    ligandContact = [i[1] for i in contactData]

    prevLigandContact = 0
    digitalChange = []
    for nextLigandContact in ligandContact:
        if nextLigandContact != prevLigandContact:
            digitalChange.append(1)
            prevLigandContact = nextLigandContact
        else:
            digitalChange.append(0)

    totalLifeTime = 0

    stablePacketALL = []
    stablePacket = []
    for digitalChange_i in digitalChange:
        if digitalChange_i == 0:
            stablePacket.append(0)
            totalLifeTime += 1
        elif digitalChange_i == 1:
            if stablePacket:
                stablePacketALL.append(stablePacket)
                stablePacket = []
        else:
            raise ValueError
    stablePacketALL.append(stablePacket)

    averageLifeTime = np.average([len(stablePack) for stablePack in stablePacketALL])
    maxLifeTime = np.max([len(stablePack) for stablePack in stablePacketALL])
    minLifeTime = np.min([len(stablePack) for stablePack in stablePacketALL])    

    return totalLifeTime, averageLifeTime, maxLifeTime, minLifeTime

def cleanContactData(contactData, budTime):
    contactDataTRIMMED = []
    for contactData_i in contactData:
        if contactData_i[0] <= budTime:
            contactDataTRIMMED.append(contactData_i)
    return contactDataTRIMMED
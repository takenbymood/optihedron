import networkx as nx
import random
import numpy as np
import nanoparticle
import run
import os

def gnmw(G):
    return(len(G.nodes()),len(G.edges()),np.sum([w['weight'] for (u,v,w) in G.edges(data=True)]))

def gnmRandomWeightedGraph(nodes,edges,weight):
    G = nx.gnm_random_graph(nodes,edges)
    totalWeight = 0
    for (u,v,w) in G.edges(data=True):
        w['weight'] = random.uniform(0,10)
        totalWeight += w['weight']
    
    normFactor = weight/totalWeight
    
    for (u,v,w) in G.edges(data=True):
        w['weight'] = w['weight']*normFactor
        
    return G

def buildNanoParticleFromNetwork(G,weight,radius=4,sig=1):
    particle = nanoparticle.NanoParticle()
    for n,w in G.nodes(data=True):
        particle.ligands.append(nanoparticle.Ligand(weight,sig,radius,w['pol'],w['azi']))
    return particle

def gArcDist(Ang1, Ang2, rad=4):
    #Ang = (PolarAng,AziAng)
    #https://math.stackexchange.com/questions/231221/great-arc-distance-between-two-points-on-a-unit-sphere
    arcDist=rad*(np.arccos((np.cos(Ang1[0])*np.cos(Ang2[0]))+((np.sin(Ang1[0])*np.sin(Ang2[0]))*(np.cos(Ang1[1]-Ang2[1])))))
    return arcDist

def createSmallWorldLigandNetwork(n,p):
    SWG = nx.watts_strogatz_graph(n,n-1,p, seed=None)
    minWeight = 100
    maxWeight = 0
    azipol = []
    for (n,w) in SWG.nodes(data=True):
        delta = 0
        azi = 0
        pol = 0
        while delta < 1.65:
            azi = random.uniform(0,6.2831)
            pol = random.uniform(0,3.141)
            minDelta = 100
            for n in azipol:
                d = gArcDist((pol,azi),(n[0],n[1]))
                if d < minDelta:
                    minDelta = d
            delta = minDelta
        azipol.append((pol,azi))
        w['azi'] = azi
        w['pol'] = pol

    for (u,v,w) in SWG.edges(data=True):
        unode = SWG.nodes(data=True)[u]
        vnode = SWG.nodes(data=True)[v]
        w['weight'] = 1.0/gArcDist((unode['pol'], unode['azi']),(vnode['pol'], vnode['azi']))
    return SWG

def smallWorldNess(G):
    gnmwG = gnmw(G)
    RG = gnmRandomWeightedGraph(gnmwG[0],gnmwG[1],gnmwG[2])
    pLengthG = nx.average_shortest_path_length(G,weight='weight')
    pLengthRG = nx.average_shortest_path_length(RG,weight='weight')
    clusteringG = nx.average_clustering(G,weight='weight')
    clusteringRG = nx.average_clustering(RG,weight='weight')
    SW = (clusteringG/clusteringRG)/(pLengthG/pLengthRG)
    return SW

def runSmallWorldSweep(n,nodes,weight):
    metrics = []
    for i in range(n):
        net = createSmallWorldLigandNetwork(nodes,0.1)
        s = smallWorldNess(net)
        particle = buildNanoParticleFromNetwork(net,weight)
        f = run.evaluateParticle(particle,"sweep_"+str(i))
        metrics.append((s,f))
    with open(os.path.join('smallworld.csv'), 'a') as file_:
        for m in metrics:
            file_.write(str(m).replace('(','').replace(')','').replace(' ','')[:-2]+"\n")


runSmallWorldSweep(100,25,8)

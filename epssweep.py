import networkx as nx
import random
import numpy as np
import nanoparticle
import run
import os


def gArcDist(Ang1, Ang2, rad=4):
    #Ang = (PolarAng,AziAng)
    #https://math.stackexchange.com/questions/231221/great-arc-distance-between-two-points-on-a-unit-sphere
    arcDist=rad*(np.arccos((np.cos(Ang1[0])*np.cos(Ang2[0]))+((np.sin(Ang1[0])*np.sin(Ang2[0]))*(np.cos(Ang1[1]-Ang2[1])))))
    return arcDist

def buildNanoParticleFromNetwork(G,weight,radius=4,sig=1):
    particle = nanoparticle.NanoParticle()
    for n,w in G.nodes(data=True):
        particle.ligands.append(nanoparticle.Ligand(weight,sig,radius,w['pol'],w['azi']))
    return particle

def createRandomLigandNetwork(n):
    SWG = nx.complete_graph(n)
    minWeight = 100
    maxWeight = 0
    azipol = []
    for (n,w) in SWG.nodes(data=True):
        delta = 0
        azi = 0
        pol = 0
        attempts = 0
        while delta < 1.65 and attempts < 100:
            azi = random.uniform(0,6.2831)
            pol = random.uniform(0,3.141)
            minDelta = 100
            for n in azipol:
                d = gArcDist((pol,azi),(n[0],n[1]))
                if d < minDelta:
                    minDelta = d
            delta = minDelta
            attempts += 1
        azipol.append((pol,azi))
        w['azi'] = azi
        w['pol'] = pol
        
    if attempts >= 100:
            return createRandomLigandNetwork(n)

    for (u,v,w) in SWG.edges(data=True):
        unode = SWG.nodes(data=True)[u]
        vnode = SWG.nodes(data=True)[v]
        w['weight'] = 1.0/gArcDist((unode['pol'], unode['azi']),(vnode['pol'], vnode['azi']))
    return SWG

def runEpsSweep(minNodes,maxNodes,totalEps):
    metrics = []
    for nodes in range(minNodes,maxNodes+1):
        weight = float(totalEps)/float(nodes)
        net = createRandomLigandNetwork(nodes)
        particle = buildNanoParticleFromNetwork(net,weight)
        f = run.evaluateParticle(particle,"sweep_"+str(nodes))
        metrics.append((nodes,f))
        with open(os.path.join('epssweep_'+str(totalEps)+'.csv'), 'a') as file_:
            file_.write(str((nodes,f)).replace('(','').replace(')','').replace(' ','')[:-2]+"\n")


runEpsSweep(20,72,450)
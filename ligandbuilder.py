import networkx as nx
import random
import numpy as np
import nanoparticle
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

def createLigandNetworkFromGraph(G):
    azipol = []
    for (n,w) in G.nodes(data=True):
        delta = 0
        azi = 0
        pol = 0
        attempts = 0
        while delta < 1.65 and attempts < 1000:
            azi = random.uniform(0,6.2831)
            pol = random.uniform(0,3.141)
            minDelta = 100
            for n in azipol:
                d = gArcDist((pol,azi),(n[0],n[1]))
                if d < minDelta:
                    minDelta = d
            delta = minDelta
            attempts += 1
        if attempts >= 100:
            break
        azipol.append((pol,azi))
        w['azi'] = azi
        w['pol'] = pol
        
    if attempts >= 100:
        return createLigandNetworkFromGraph(G)
    else:
        for (u,v,w) in G.edges(data=True):
            unode = G.nodes(data=True)[u]
            vnode = G.nodes(data=True)[v]
            w['weight'] = 1.0/gArcDist((unode['pol'], unode['azi']),(vnode['pol'], vnode['azi']))
    return G

def createRandomLigandNetwork(n):
    G = nx.complete_graph(n)
    G = createLigandNetworkFromGraph(G)
    return G

def createSmallWorldLigandNetwork(n,p):
    SWG = nx.watts_strogatz_graph(n,n-1,p, seed=None)
    G = createLigandNetworkFromGraph(SWG)
    return SWG
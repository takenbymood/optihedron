import networkx as nx
import random
import numpy as np
import nanoparticle
import run
import os
import ligandbuilder as lg
import smallworld as sw

def runParticle(particle,simName):
    f = run.evaluateParticle(particle,simName)


def runEpsSweep(minNodes,maxNodes,totalEps):
    metrics = []
    for nodes in range(minNodes,maxNodes+1):
        weight = float(totalEps)/float(nodes)
        net = lg.createRandomLigandNetwork(nodes)
        particle = lg.buildNanoParticleFromNetwork(net,weight)
        f = runParticle(particle,"sweep_"+str(nodes))
        metrics.append((nodes,f))
        with open(os.path.join('epssweep_'+str(totalEps)+'.csv'), 'a') as file_:
            file_.write(str((nodes,f)).replace('(','').replace(')','').replace(' ','')[:-2]+"\n")

def runSmallWorldSweep(n,nodes,weight):
    metrics = []
    for i in range(n):
        net = lg.createSmallWorldLigandNetwork(nodes,0.1)
        s = sw.smallWorldNess(net)
        particle = lg.buildNanoParticleFromNetwork(net,weight)
        f = runParticle(particle,"sweep_"+str(i))
        metrics.append((s,f))
    with open(os.path.join('smallworld_'+str(nodes)+'_'+str(weight)+'.csv'), 'a') as file_:
        for m in metrics:
            file_.write(str(m).replace('(','').replace(')','').replace(' ','')[:-2]+"\n")


runEpsSweep(30,30,150)
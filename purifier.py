from db import databaseconnection as dbc
import sqlalchemy
from sqlalchemy import Table, create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import *
from sqlalchemy.orm import *
import random
import os
import time
import nanoparticle
import networkx as nx

import numpy as np
from ga import networkedgeneticalgorithm as nga
from nanoparticle import Ligand, NanoParticle

from tools import analysistools as atools

import colldb
from colldb import Particle, Instance, indToParticle, sessToInst

import argparse

import sys


parser = argparse.ArgumentParser(description='')

parser.add_argument('-out','--outdir', default="", type=str, 
                    help='output directory for the generated csv files')
parser.add_argument('-db','--database', default='', type=str, 
                    help='input database, must be in sqlite format')

args = parser.parse_args()

dbPath = args.database

outdir = args.outdir

Base = declarative_base()
engine = sqlalchemy.create_engine('sqlite:///{}'.format(dbPath))
Base.metadata.create_all(bind=engine)
dbSession = sessionmaker(bind=engine)
dbSession = dbSession()
budFilePath = os.path.join(outdir,'opti-bud.csv')
ffFilePath = os.path.join(outdir,'opti-ff.csv')

netFilePaths = []

minPrune = 0.0
maxPrune = 1.0
pruneStep = 0.1

pruneSteps = [x for x in atools.frange(minPrune,maxPrune,pruneStep)]

for pruning in pruneSteps:
    netFilePaths.append(os.path.join(outdir,'opti-net'+str(pruning)+'.csv'))


for n in netFilePaths:
    with open(n, 'w') as netFile:
        netWriter = atools.UnicodeWriter(netFile)
        netWriter.writerows([[
            "Ligand Number",
            "Mean Affinity",
            "Fitness",
            "Mean Budding Time",
            "Budding Rate",
            "Density",
            "Max Diameter",
            "Mean Diameter",
            "Min Diameter",
            "Min Radius",
            "Mean Radius",
            "Max Radius",
            "Max Average Shortest Path",
            "Mean Average Shortest Path",
            "Min Average Shortest Path",
            "Subgraph Number",
            "Estrada Coefficient",
            "Pruning",
            "Max SmallWorld",
            "Mean SmallWorld",
            "Min SmallWorld"
        ]])

with open(ffFilePath, 'w') as ffFile:
    ffWriter = atools.UnicodeWriter(ffFile)
    ffWriter.writerows([[
        "Ligand Number",
        "Mean Affinity",
        "Fitness",
        "Mean Budding Time",
        "Budding Rate",
        "Patchiness",
        "Lininess",
        "Spottiness",
        "Walks",
        "End To End Distance",
        "Furthest Ligands"
        ]])


    buddingRates = {}
    buddingTimes = {}
    sizes = {}
    budded = {}

    budIndividuals = 0
    noBudIndividuals = 0
    printed = []
    for i in dbSession.query(Particle).yield_per(100):

        nL = int(float(i.nligands))
        aE = int(np.round(float(i.avgEps)))
        if not (nL,aE) in printed:
            printed.append((nL,aE))
            print((nL,aE))
        if not nL in buddingRates:
            buddingRates[nL] = {}
            buddingTimes[nL] = {}
            sizes[nL] = {}
            budded[nL] = {}
        if not aE in buddingRates[nL]:
            buddingRates[nL][aE] = 0
            buddingTimes[nL][aE] = 0
            sizes[nL][aE] = 0
            budded[nL][aE] = 0

        buddingRates[nL][aE] += i.budPerc
        buddingTimes[nL][aE] += i.budTime
        sizes[nL][aE] += 1
        if i.budTime > 0.0:
        	budded[nL][aE]+=1

        # if i.walks != None:
        #     walkStrs = ""
        #     for walk in i.walks:
        #         walkStr = "{"
        #         stepCount = 0
        #         for s in walk.trajectory:
        #             if len(s[1])>0 or len(s[2])>0:
        #                 walkStr+="s" + str(stepCount) + '['
        #                 for c in s[1]:
        #                     walkStr += 'c' + str(c)
        #                     walkStr += ';'
        #                 for d in s[2]:
        #                     walkStr += 'd' + str(d)
        #                     walkStr += ';'
        #                 if walkStr[-1] == ';':
        #                     walkStr = walkStr[:-1]
        #                 walkStr += ']'
        #             stepCount += 1
        #         walkStr+="}"
        #         walkStrs += walkStr

        lNum = 0
        dists = {}
        for l in i.phenome.particle.ligands:
            if l.eps > 0:
                dists[lNum] = {}
                mNum = 0
                for m in i.phenome.particle.ligands:
                    if m.eps > 0:
                        dists[lNum][mNum] = np.sqrt(np.sum([x*x for x in np.subtract(atools.sphPol2Crt((l.rad,l.polAng,l.aziAng)),atools.sphPol2Crt((m.rad,m.polAng,m.aziAng)))]))
                    mNum += 1 
            lNum+=1
        maxDist = 0
        maxPair = ()
        for k,v in dists.iteritems():
            for k2,v2 in v.iteritems():
                if v2 > maxDist:
                    maxDist = v2
                    maxPair = (k,k2)



        ffWriter.writerows([[
            str(i.nligands),
            str(i.avgEps),
            str(i.fitness),
            str(i.budTime),
            str(i.budPerc),
            str(i.patchiness),
            str(i.lininess),
            str(i.spottiness),
            '',
            str(maxDist),
            '('+str(maxPair[0])+';'+str(maxPair[1])+')'
            ]])
        
        c = 0
        for n in netFilePaths:
            with open(n, 'a') as netFile:
                netWriter = atools.UnicodeWriter(netFile)
                pN = atools.pruneNetwork(i.network,pruneSteps[c])

                density = nx.density(pN)
                graphs = list(nx.connected_component_subgraphs(pN))
                dS = []
                rS = []
                sPS = []
                SWs = []

                for g in graphs:
                    d = nx.diameter(g)
                    r = nx.radius(g)
                    sp = nx.average_shortest_path_length(g)
                    dS.append(d)
                    rS.append(r)
                    sPS.append(sp)
                    SWs.append(atools.smallWorldNess(g))

                maxDiameter = np.max(dS)
                avgDiameter = np.mean(dS)
                minDiameter = np.min(dS)
                minRadius = np.min(rS)
                avgRadius = np.mean(rS)
                maxRadius = np.max(rS)
                maxASp = np.max(sPS)
                avgASp = np.mean(sPS)
                minASp = np.min(sPS)
                maxSW = np.max(SWs)
                avgSW = np.mean(SWs)
                minSW = np.min(SWs)

                subgraphs = len(graphs)

                estrada = nx.estrada_index(pN)

                netWriter.writerows([[
                    str(i.nligands),
                    str(i.avgEps),
                    str(i.fitness),
                    str(i.budTime),
                    str(i.budPerc),
                    str(density),
                    str(maxDiameter),
                    str(avgDiameter),
                    str(minDiameter),
                    str(minRadius),
                    str(avgRadius),
                    str(maxRadius),
                    str(maxASp),
                    str(avgASp),
                    str(minASp),
                    str(subgraphs),
                    str(estrada),
                    str(pruneSteps[c]),
                    str(maxSW),
                    str(avgSW),
                    str(minSW)
                    ]])
            c+=1

        if i.budTime != -1:
            budIndividuals += 1
        else:
            noBudIndividuals += 1
       
    for k,v in buddingRates.iteritems():
        for k2, v2 in v.iteritems():
            sn = float(sizes[k][k2])
            bn = float(budded[k][k2])
            avgBudding = 0.0
            avgBudTime = -1.0
            if sn > 0.0:
                avgBudding = float(v2)/sn
            if bn > 0.0:
            	avgBudTime = float(buddingTimes[k][k2])/bn

            buddingRates[k][k2] = avgBudding
            buddingTimes[k][k2] = avgBudTime

with open(budFilePath, 'w') as budFile:
    budWriter = atools.UnicodeWriter(budFile)
    for k,v in buddingRates.iteritems():
        for k2, v2 in v.iteritems():
            budWriter.writerows([[
                str(k),
                str(k2),
                str(budded[k][k2]),
                str(sizes[k][k2]),
                str(v2),
                str(buddingTimes[k][k2])
                ]])

    print("Budded: " + str(budIndividuals))
    print("Did Not Bud: " + str(noBudIndividuals))
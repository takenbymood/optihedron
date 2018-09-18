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

import sys

dbPath = '/Users/joel/Projects/optidb/sep.db'

Base = declarative_base()
engine = sqlalchemy.create_engine('sqlite:///{}'.format(dbPath))
Base.metadata.create_all(bind=engine)
dbSession = sessionmaker(bind=engine)
dbSession = dbSession()
budFilePath = '/Users/joel/Projects/optidb/opti-bud.csv'
ffFilePath = '/Users/joel/Projects/optidb/opti-ff.csv'
netFilePath = '/Users/joel/Projects/optidb/opti-net.csv'
with open(ffFilePath, 'w') as ffFile, open(netFilePath, 'w') as netFile:
    ffWriter = atools.UnicodeWriter(ffFile)
    netWriter = atools.UnicodeWriter(netFile)

    ffWriter.writerows([[
        "Ligand Number",
        "Average Affinity",
        "Fitness",
        "Average Budding Time",
        "Budding Rate",
        "Patchiness",
        "Lininess",
        "Spottiness"
        ]])

    netWriter.writerows([[
        "Ligand Number",
        "Average Affinity",
        "Fitness",
        "Average Budding Time",
        "Budding Rate",
        "Density",
        "Max Diameter",
        "Mean Diameter",
        "Min Radius",
        "Average Radius",
        "Subgraph Number",
        "Estrada Coefficient",
        "Rich Club Coefficients"
        ]])

    buddingRates = {}
    buddingTimes = {}
    sizes = {}
    budded = {}

    budIndividuals = 0
    noBudIndividuals = 0
    for i in dbSession.query(Particle).yield_per(100).limit(5000):

        nL = int(float(i.nligands))
        aE = int(np.round(float(i.avgEps)))
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

        ffWriter.writerows([[
            str(i.nligands),
            str(i.avgEps),
            str(i.fitness),
            str(i.budTime),
            str(i.budPerc),
            str(i.patchiness),
            str(i.lininess),
            str(i.spottiness)
            ]])
        
        pN = atools.pruneNetwork(i.network,0.5)

        density = nx.density(pN)
        graphs = list(nx.connected_component_subgraphs(pN))
        dS = []
        rS = []
        for g in graphs:
            d = nx.diameter(g)
            r = nx.radius(g)
            dS.append(d)
            rS.append(r)
        maxDiameter = np.max(dS)
        avgDiameter = np.mean(dS)
        minRadius = np.min(rS)
        avgRadius = np.mean(rS)
        subgraphs = len(graphs)

        estrada = nx.estrada_index(pN)

        richclub = nx.rich_club_coefficient(pN,normalized=False)

        netWriter.writerows([[str(i.nligands),
            str(i.avgEps),
            str(i.fitness),
            str(i.budTime),
            str(i.budPerc),
            str(density),
            str(maxDiameter),
            str(avgDiameter),
            str(minRadius),
            str(avgRadius),
            str(subgraphs),
            str(estrada),
            str(richclub)
            ]])

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
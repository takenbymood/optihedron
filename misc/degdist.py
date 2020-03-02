from tools import analysistools as atools
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import random
from scipy.stats import sem
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.spatial import distance_matrix
import networkx as nx

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

import pickle
import nanoparticle

from operator import add

import numpy as np
from ga import networkedgeneticalgorithm as nga
from nanoparticle import Ligand, NanoParticle

import colldb
from colldb import Particle, Instance, indToParticle, sessToInst

import argparse

import sys
import seaborn as sns
import collections

cR = [[20, 13], [20, 14], [21, 12], [21, 13], [21, 14], [22, 12], [22, 13], [22, 14], [23, 10], [23, 11], [23, 12], [23, 13], [24, 10], [24, 11], [24, 12], [24, 13], [25, 9], [25, 10], [25, 11], [25, 12], [26, 9], [26, 10], [26, 11], [27, 9], [27, 10], [28, 9], [28, 10], [29, 8], [29, 9], [30, 8], [30, 9], [31, 7], [31, 8], [32, 7], [32, 8], [33, 7], [33, 8], [34, 7], [35, 6], [35, 7], [36, 6], [37, 6], [38, 6], [39, 6], [40, 6], [42, 5], [43, 5], [44, 5], [45, 5], [46, 5], [50, 4], [51, 4], [52, 4], [53, 4], [54, 4]]


dbPath = "/Users/joelforster/Projects/optidb/sep2.db"
Base = declarative_base()
engine = sqlalchemy.create_engine('sqlite:///{}'.format(dbPath))
Base.metadata.create_all(bind=engine)
dbSession = sessionmaker(bind=engine)
dbSession = dbSession()

dbSession.expire_all()
dbSession.flush()

degDists = {}

printed = []
s = 0
for i in dbSession.query(Particle).yield_per(100):
    nL = int(float(i.nligands))
    aE = int(np.round(float(i.avgEps)))

    if not [nL,aE] in cR:
        continue
    if not nL in degDists:
        degDists[nL] = {}
    if not aE in degDists[nL]:
        degDists[nL][aE] = {}
        degDists[nL][aE]['b'] = {}
        degDists[nL][aE]['n'] = {}
        degDists[nL][aE]['b']['num'] = 0
        degDists[nL][aE]['n']['num'] = 0
    if not (nL,aE) in printed:
        printed.append((nL,aE))
        print((nL,aE))
    G = atools.pruneNetwork(i.network,0.3)
    degrees=sorted([d for n, d in G.degree()], reverse=True)
    degreeCount = collections.Counter(degrees)
    deg, cnt = zip(*degreeCount.items())
    degPad = []
    cntPad = []
    for n in range(nL):
        if not n in deg:
            degPad.append(n)
            cntPad.append(0)
        else:
            degPad.append(n)
            cntPad.append(cnt[deg.index(n)])
    if i.budTime > 0:
        if not 'deg' in degDists[nL][aE]['b']:
            degDists[nL][aE]['b']['deg'] = degPad
            degDists[nL][aE]['b']['cnt'] = cntPad
        degDists[nL][aE]['b']['num'] += 1
        degDists[nL][aE]['b']['cnt'] = list(map(add,degDists[nL][aE]['b']['cnt'], cntPad) )
    else:
        if not 'deg' in degDists[nL][aE]['n']:
            degDists[nL][aE]['n']['deg'] = degPad
            degDists[nL][aE]['n']['cnt'] = cntPad
        degDists[nL][aE]['n']['num'] += 1
        degDists[nL][aE]['n']['cnt'] = list(map(add,degDists[nL][aE]['n']['cnt'], cntPad) )
    s+=1

for nL,v1 in degDists.iteritems():
    for aE,v in v1.iteritems():
        if v['b']['num'] > 0:
            v['b']['avg'] = [float(d)/float(v['b']['num']) for d in v['b']['cnt']]
        if v['n']['num'] > 0:
            v['n']['avg'] = [float(d)/float(v['n']['num']) for d in v['n']['cnt']]

sns.set_color_codes("pastel")

for nL,v1 in degDists.iteritems():
    for aE,v in v1.iteritems():
        if ('avg' in v['b']) and ('avg' in v['n']):
            bdf = pd.DataFrame()
            bdf['nlg'] = [nL for d in v['b']['avg']]
            bdf['aE'] = [aE for d in v['b']['avg']]
            bdf['deg'] = v['b']['deg']
            bdf['avg'] = v['b']['avg']
            bdf['cat'] = ['budding' for d in v['b']['avg']]

            ndf = pd.DataFrame()
            ndf['nlg'] = [nL for d in v['n']['avg']]
            ndf['aE'] = [aE for d in v['n']['avg']]
            ndf['deg'] = v['n']['deg']
            ndf['avg'] = v['n']['avg']
            ndf['cat'] = ['non budding' for d in v['n']['avg']]

            df = pd.DataFrame()
            df = bdf.append(ndf)

            degDists[nL][aE]['df'] = df

            fig, ax = plt.subplots()

            sns.barplot(x='deg', y ='avg',hue='cat',data=df)
            plt.xlim(-1, max(df.query('avg > 0.0')['deg'].values)+1)

            plt.title("Mean Degree Distribution ("+str(nL)+","+str(aE)+")")
            plt.ylabel("Count")
            plt.xlabel("Node Degree")
            plt.legend()
            plt.savefig("/Users/joelforster/Projects/optidb/clus/degs/"+"("+str(nL)+","+str(aE)+")")

bigDF = pd.DataFrame()

for nL,v1 in degDists.iteritems():
    for aE,v in v1.iteritems():
        if 'df' in v:
            bigDF = bigDF.append(v['df'])

bigDF.to_csv("/Users/joelforster/Projects/optidb/clus/degs/df.csv",index=False)

bigDF = None

with open('/Users/joelforster/Projects/optidb/clus/degdist.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(degDists, f, pickle.HIGHEST_PROTOCOL)

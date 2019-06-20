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


cR = [[20, 13], [20, 14], [21, 12], [21, 13], [21, 14], [22, 12], [22, 13], [22, 14], [23, 10], [23, 11], [23, 12], [23, 13], [24, 10], [24, 11], [24, 12], [24, 13], [25, 9], [25, 10], [25, 11], [25, 12], [26, 9], [26, 10], [26, 11], [27, 9], [27, 10], [28, 9], [28, 10], [29, 8], [29, 9], [30, 8], [30, 9], [31, 7], [31, 8], [32, 7], [32, 8], [33, 7], [33, 8], [34, 7], [35, 6], [35, 7], [36, 6], [37, 6], [38, 6], [39, 6], [40, 6], [42, 5], [43, 5], [44, 5], [45, 5], [46, 5], [50, 4], [51, 4], [52, 4], [53, 4], [54, 4]]


args = parser.parse_args()

dbPath = args.database

outdir = args.outdir

Base = declarative_base()
engine = sqlalchemy.create_engine('sqlite:///{}'.format(dbPath))
Base.metadata.create_all(bind=engine)
dbSession = sessionmaker(bind=engine)
dbSession = dbSession()
spectrumFilePath = os.path.join(outdir,'opti-spectrum.csv')
spectra = []

with open(spectrumFilePath, 'w') as specFile:
	specWriter = atools.UnicodeWriter(specFile)
	specWriter.writerows([['id','e','n','fitness','budtime','budfrac','mean_dist','max_dist','min_dist','spectrum']])

for i in dbSession.query(Particle).yield_per(100):
	nL = int(float(i.nligands))
	aE = int(np.round(float(i.avgEps)))

	if not [nL,aE] in cR:
		continue
	nodelist = []
	for n in i.network.nodes(data=True):
		if n[1]['weight'] != 0:
			nodelist.append(n[0])
	adj = nx.to_numpy_matrix(i.network,nodelist=nodelist)
	spectrum = np.linalg.eigvals(adj)
	mask = np.ones(adj.shape, dtype=bool)
	np.fill_diagonal(mask, 0)
	meanDist = adj[mask].mean()
	minDist = adj[mask].min()
	maxDist = adj[mask].max()
	spectra.append([str(int(i.pID)),str(float(i.avgEps)),str(float(i.nligands)),str(float(i.fitness)),str(float(i.budTime)),str(float(i.budPerc)),str(float(meanDist)),str(float(maxDist)),str(float(minDist)),str(list(spectrum))])
	if len(spectra) == 100:
		with open(spectrumFilePath, 'a') as specFile:
			specWriter = atools.UnicodeWriter(specFile)
			specWriter.writerows(spectra)
			spectra = []
from db import databaseconnection as dbc
import sqlalchemy
from sqlalchemy import Table, create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref
import random
import os
import time

from membranesimulation import MembraneSimulation

from tools import vectools

import numpy as np
from ga import networkedgeneticalgorithm as nga
from nanoparticle import Ligand, NanoParticle

import colldb
from colldb import Particle, Instance, indToParticle, sessToInst

import argparse

import sys


parser = argparse.ArgumentParser(description='')

parser.add_argument('-out','--outdir', default="", type=str, 
                    help='output directory for the generated csv files')
parser.add_argument('-db','--database', default='', type=str, 
                    help='input database, must be in sqlite format')

parser.add_argument('-nL','--ligands', default=-1, type=int, 
                    help='number of ligands')
parser.add_argument('-aE','--epsilon', default=-1, type=int, 
                    help='particle affinity')

args = parser.parse_args()

dbPath = args.database

outdir = args.outdir

ligandNum = args.ligands
avgAffinity = args.epsilon

PDIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATEDIR = os.path.join(PDIR,'mem/template')
TEMPLATEDATAPATH = os.path.join(TEMPLATEDIR,'data.template')
TEMPLATEINPUTPATH = os.path.join(TEMPLATEDIR,'in.template')

Base = declarative_base()
engine = sqlalchemy.create_engine('sqlite:///{}'.format(dbPath))
Base.metadata.create_all(bind=engine)
dbSession = sessionmaker(bind=engine)
dbSession = dbSession()

cR = [[20, 13], [20, 14], [21, 12], [21, 13], [21, 14], [22, 12], [22, 13], [22, 14], [23, 10], [23, 11], [23, 12], [23, 13], [24, 10], [24, 11], [24, 12], [24, 13], [25, 9], [25, 10], [25, 11], [25, 12], [26, 9], [26, 10], [26, 11], [27, 9], [27, 10], [28, 9], [28, 10], [29, 8], [29, 9], [30, 8], [30, 9], [31, 7], [31, 8], [32, 7], [32, 8], [33, 7], [33, 8], [34, 7], [35, 6], [35, 7], [36, 6], [37, 6], [38, 6], [39, 6], [40, 6], [42, 5], [43, 5], [44, 5], [45, 5], [46, 5], [50, 4], [51, 4], [52, 4], [53, 4], [54, 4]]
# .filter(int(float(Particle.avg_eps))==avgAffinity)

for i in dbSession.query(Particle).filter_by(nligands=ligandNum).filter_by(avgEps=avgAffinity).yield_per(100):
        nL = int(float(i.nligands))
        aE = int(np.round(float(i.avgEps)))

        sim = MembraneSimulation(
            str(i.pID)+"_"+str(int(i.budTime))+"_"+str(int(i.budPerc)),
            i.phenome.particle,
            25000,
            0.01,        
            outdir,
            outdir,
            TEMPLATEDATAPATH,
            TEMPLATEINPUTPATH,
            rAxis=vectools.randomUnitVector(),
            rAmount=random.uniform(0.3141,3.141)
        )
        sim.saveFiles()
        

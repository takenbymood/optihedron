from db import databaseconnection as dbc
from sqlalchemy import Table, create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy import *
from sqlalchemy.orm import *
import random
import os
import time
import nanoparticle

import numpy as np
from ga import networkedgeneticalgorithm as nga
from nanoparticle import Ligand, NanoParticle

from tools import analysistools as atools

dbParentPath = '/Users/joelforster/Projects/optidb/sep'
dbPaths = os.listdir(dbParentPath)

import sys

def progressbar(it, prefix="", size=60):
    count = len(it)
    def _show(_i):
        x = int(size*_i/count)
        sys.stdout.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), _i, count))
        sys.stdout.flush()

    _show(0)
    for i, item in enumerate(it):
        yield item
        _show(i+1)
    sys.stdout.write("\n")
    sys.stdout.flush()

Base = declarative_base()

class Instance(Base):
    __tablename__ ='instances'
    
    pID = Column('id', Integer, primary_key=True)
    timestamp = Column('timestamp', String)
    originalDB = Column('original_db', String)
    arguments = Column('arguments', String)
    
    particles = relationship('Particle',back_populates='instance')
    
class Particle(Base):
    __tablename__ = 'particles'

    pID = Column('id', Integer, primary_key=True)
    
    inst_id = Column(Integer, ForeignKey('instances.id'))
    instance = relationship('Instance',uselist=False)
    
    gen = Column('gen', Integer)
    fitness = Column('fitness', Numeric)
    genome = Column('genome',String)
    phenome = Column('phenome', PickleType)
    avgEps = Column('avg_eps',Numeric)
    budPerc = Column('budding_rate',Numeric)
    budTime = Column('bud_time',Numeric)
    network = Column('network', PickleType)
    patchiness = Column('patchiness', Numeric)
    lininess = Column('lininess', Numeric)
    spottiness = Column('spottiness', Numeric)
    nligands = Column('nligands', Integer)


def sessToInst(session):
    inst = Instance()
    inst.timestamp = session.timestamp
    inst.arguments = session.arguments
    return inst

def indToParticle(ind):
    part = Particle()
    part.gen = ind.gen.genNumber
    part.fitness = ind.fitness
    part.genome = ind.genome
    part.phenome = ind.phenomePickle
    part.budPerc = ind.budPerc
    part.budTime = ind.budTime
    part.network = atools.buildLigandNetwork(ind.phenomePickle.particle.ligands)
    ff = atools.formFactor(ind.phenomePickle.particle.ligands)
    part.patchiness = ff[0]
    part.lininess = ff[1]
    part.spottiness = ff[2]
    lNum = 0.0
    lEps = 0.0
    for l in ind.phenomePickle.particle.ligands:
        if l.eps>0.0:
            lNum+=1
            lEps+=l.eps
    lAvg = float(lEps)/float(lNum)
    part.nligands = lNum
    part.avgEps = lAvg
    return part


cDB = '/Users/joelforster/Projects/optidb/sep.db'
engine = create_engine('sqlite:///{}'.format(cDB))
Base.metadata.create_all(bind=engine)
dbSession = sessionmaker(bind=engine)
dbs= dbSession()

for p in dbPaths:
    if p.split('.')[-1] != "db" or p.split('-')[0]=='eps10' :
        continue
    path = os.path.join(dbParentPath,p)
    print(path)
    c = dbc.DatabaseConnection(path)
    s = c.getLastSession()
    print(s)
    time.sleep(1)
    instance = sessToInst(s)
    instance.originalDB = p
    for i in progressbar(s.getIndividualsList()):
        instance.particles.append(indToParticle(i))
    dbs.add(instance)
    dbs.commit()
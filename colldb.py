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


import sys

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

    walks = relationship('Walk',back_populates='particle')
    
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

class Walk(Base):
    __tablename__ = 'walks'
    pID = Column('id', Integer, primary_key=True)
    
    particle_id = Column(Integer, ForeignKey('particles.id'))
    particle = relationship('Particle',uselist=False)

    trajectory = Column('trajectory', PickleType)


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

    for sim in ind.sims:
        if len(sim.data) >0 and 'walk' in sim.data[-1]:
            walk = Walk()
            walk.trajectory = []
            for s in sim.data:
                walk.trajectory.append(s['walk'])
            part.walks.append(walk)

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
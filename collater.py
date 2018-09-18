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

import colldb
from colldb import Particle, Instance, indToParticle, sessToInst

import sys

Base = declarative_base()


def main():

    dbParentPath = '/Users/joelforster/Projects/optidb/sep'
    dbPaths = os.listdir(dbParentPath)
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
        for i in atools.progressbar(s.getIndividualsList()):
            instance.particles.append(indToParticle(i))
        dbs.add(instance)
        dbs.commit()

if __name__ == "__main__":
    main()
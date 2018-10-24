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

import argparse
import sys

parser = argparse.ArgumentParser(description='')

parser.add_argument('-dir','--dir', default=None, type=str, 
                    help='base directory of all of the databases, they must by sqlite files and have the .db extension to be included')
parser.add_argument('-db','--database', default='collated.db', type=str, 
                    help='output database, will be in sqlite format')





def main():

    args = parser.parse_args()

    dbParentPath = args.dir
    dbPaths = os.listdir(dbParentPath)
    cDB = args.database

    if dbParentPath == None:
        print('no database collection specified, use the -dir argument')
        return

    engine = create_engine('sqlite:///{}'.format(cDB))
    colldb.Base.metadata.create_all(bind=engine)
    dbSession = sessionmaker(bind=engine)
    dbs= dbSession()

    for p in dbPaths:
        if p.split('.')[-1] != "db":
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
from sqlalchemy import create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref
from sqlalchemy.types import TypeDecorator

import time
import json

Base = declarative_base()

class Sessions(Base):
	__tablename__ ='sessions'

	sessionId = Column('sessionId', String, primary_key=True)	
	sessionMetrics = relationship('Metrics', backref=backref('sessions', uselist=False))
	sessionGenealogy = relationship('Genealogy', backref=backref('sessions', uselist=False))
	sessionIndividuals = relationship('Individual', backref='sessions')

	def __init__(self, sessionId):
		self.sessionId = sessionId

class Metrics(Base):
	__tablename__ = 'metrics'
	
	sessionId = Column(String, ForeignKey('sessions.sessionId'), primary_key=True)
	metricsPickle = Column('metricsPickle', PickleType)

	def __init__(self, sessionId, metrics):
		self.sessionId = sessionId
		self.metricsPickle = metrics		

class Genealogy(Base):
	__tablename__ = 'genealogy'

	sessionId = Column(String, ForeignKey('sessions.sessionId'), primary_key=True)
	treePickle = Column('treePickle', PickleType)
	historyPickle = Column('historyPickle', PickleType)

	def __init__(self, sessionId, tree, history):
		self.sessionId = sessionId	
		self.treePickle = tree
		self.historyPickle = history

class Individual(Base):
	__tablename__ = 'individuals'

	sessionId = Column(String, ForeignKey('sessions.sessionId'), primary_key=True)
	gen = Column('gen', Integer, primary_key=True)
	ind = Column('ind', Integer, primary_key=True)
	fitness = Column('fitness', Numeric)
	genomePickle = Column('genomePickle', PickleType)
	phenomePickle = Column('phenomePickle', PickleType)

	def __init__(self, sessionId, gen, ind, fitness, genome, phenome):
		self.sessionId = sessionId
		self.gen = gen
		self.ind = ind
		self.fitness = fitness
		self.genomePickle = genome
		self.phenomePickle = phenome

class DatabaseConnection:

	def __init__(self, dbfile):		
		engine = create_engine('sqlite:///{}'.format(dbfile))						
		Base.metadata.create_all(bind=engine)
		dbSession = sessionmaker(bind=engine)

		self.dbSession = dbSession()

		self.gaSessionId = time.strftime("%Y-%m-%d %H:%M:%S")
		self.gaSession = Sessions(self.gaSessionId)

		self.dbSession.add(self.gaSession)
		self.dbSession.commit()		

	def saveMetrics(self, metrics):				
		self.dbSession.add(Metrics(self.gaSessionId, metrics))					
		
	def saveGenealogy(self, tree, history):
		self.dbSession.add(Genealogy(self.gaSessionId, tree, history))		

	def saveIndividual(self, gen, ind, fitness, genome, phenome):
		self.dbSession.add(Individual(self.gaSessionId, gen, ind, fitness, genome, phenome))

	def commit(self):
		self.dbSession.commit()

	def close(self):
		self.dbSession.close()		

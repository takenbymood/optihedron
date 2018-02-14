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
	sessionFamily = relationship('Family', backref='sessions')

	def __init__(self, sessionId):
		self.sessionId = sessionId

class Metrics(Base):
	__tablename__ = 'metricsaa'
	
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

class Family(Base):
	__tablename__ = 'family'

	sessionId = Column(String, ForeignKey('sessions.sessionId'), primary_key=True)
	gen = Column('gen', Integer)
	fitness = Column('fitness', Numeric)
	genomePickle = Column('genomePickle', PickleType)
	phenomePickle = Column('phenomePickle', PickleType)

	def __init__(self, sessionId, gen, fitness, genome, phenome):
		self.sessionId = sessionId
		self.fitness = fitness
		self.genomePickle = genome
		self.phenomePickle = phenome

class DatabaseConnection:

	def __init__(self, dbfile):		
		engine = create_engine('sqlite:///{}'.format(dbfile), echo=True)						
		Base.metadata.create_all(bind=engine)
		dbSession = sessionmaker(bind=engine)

		self.dbSession = dbSession()

		self.gaSessionId = time.strftime("%Y-%m-%d %H:%M:%S")
		self.gaSession = Sessions(self.gaSessionId)

		self.dbSession.add(self.gaSession)
		self.dbSession.commit()		

	def saveMetrics(self, metrics):				
		self.dbSession.add(Metrics(self.gaSessionId, metrics))			
		self.dbSession.commit()				
		
	def saveGenealogy(self, tree, history):
		self.dbSession.add(Genealogy(self.gaSessionId, tree, history))
		self.dbSession.commit()

	def saveToFamily(self, gen, fitness, genome, phenome):
		self.dbSession.add(Family(self.gaSessionId, gen, fitness, genome, phenome))
		self.dbSession.commit()

	def closeConn(self):
		#self.dbSession.close
		#save the ga session name to a file

from sqlalchemy import create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref

import time

Base = declarative_base()

class Sessions(Base):
	__tablename__ ='sessions'

	sessionId = Column('sessionId', String, primary_key=True)	
	sessionMetrics = relationship('Metrics', uselist=False, backref=backref('sessions'))
	sessionGenealogy = relationship('Genealogy', uselist=False, backref=backref('sessions'))
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
		
	def saveSession(self):
		self.gaSessionId = time.strftime("%Y-%m-%d %H:%M:%S")
		self.gaSession = Sessions(self.gaSessionId)
		self.dbSession.add(self.gaSession)		

	def saveMetrics(self, metrics):				
		self.dbSession.add(Metrics(self.gaSessionId, metrics))					
		
	def saveGenealogy(self, tree, history):
		self.dbSession.add(Genealogy(self.gaSessionId, tree, history))		

	def saveIndividual(self, gen, ind, fitness, genome, phenome):
		self.dbSession.add(Individual(self.gaSessionId, gen, ind, fitness, genome, phenome))

	def whatSessions(self):
		gaSessions = self.dbSession.query(Sessions).all()	
		return [gaSession.sessionId for gaSession in gaSessions]	

	def loadSession(self, sessionId):
		gaSession = self.dbSession.query(Sessions).filter(Sessions.sessionId == sessionId).first()
		data = {}

		data['metrics'] = gaSession.sessionMetrics.metricsPickle

		data['genealogy'] = {}
		data['genealogy']['tree'] = gaSession.sessionGenealogy.treePickle
		data['genealogy']['history'] = gaSession.sessionGenealogy.historyPickle

		data['individuals'] = []
		for sessionIndividual in gaSession.sessionIndividuals:
			individual = {}
			individual['gen'] = sessionIndividual.gen
			individual['ind'] = sessionIndividual.ind
			individual['fitness'] = sessionIndividual.fitness
			individual['genome'] = sessionIndividual.genomePickle
			individual['phenome'] = sessionIndividual.phenomePickle

			data['individuals'].append(individual)

		return data

	def commit(self):
		self.dbSession.commit()

	def close(self):
		self.dbSession.close()		

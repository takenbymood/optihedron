from sqlalchemy import create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref

import time

Base = declarative_base()

class Sessions(Base):
	__tablename__ ='sessions'

	sessionId = Column('id', Integer, primary_key=True)
	sessionTimeStamp = Column('timestamp', String)	
	sessionMetrics = relationship('Metrics',uselist=False, back_populates='session')
	sessionGenealogy = relationship('Genealogy', uselist=False, back_populates='session')
	sessionIndividuals = relationship('Individual', back_populates='session')

	def __init__(self, sessionTimeStamp):
		self.sessionTimeStamp = sessionTimeStamp

class Metrics(Base):
	__tablename__ = 'metrics'
	
	pID = Column('id', Integer, primary_key=True)
	sessionId = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Sessions',uselist=False, back_populates='sessionMetrics')
	metricsPickle = Column('metricsPickle', PickleType)

	def __init__(self, sessionId, metrics):
		self.metricsPickle = metrics		

class Genealogy(Base):
	__tablename__ = 'genealogy'

	pID = Column('id', Integer, primary_key=True)
	sessionId = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Sessions',uselist=False, back_populates='sessionGenealogy')
	treePickle = Column('treePickle', PickleType)
	historyPickle = Column('historyPickle', PickleType)

	def __init__(self, sessionId, tree, history):
		self.treePickle = tree
		self.historyPickle = history

class Individual(Base):
	__tablename__ = 'individuals'

	pID = Column('id', Integer, primary_key=True)
	sessionId = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Sessions',uselist=False, back_populates='sessionIndividuals') 
	gen = Column('gen', Integer)
	ind = Column('ind', Integer)
	fitness = Column('fitness', Numeric)
	genomePickle = Column('genomePickle', PickleType)
	phenomePickle = Column('phenomePickle', PickleType)

	def __init__(self, sessionId, gen, ind, fitness, genome, phenome):
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
		self.gaSessionTimeStamp = time.strftime("%Y-%m-%d %H:%M:%S")
		self.gaSession = Sessions(self.gaSessionTimeStamp)
		self.gaSessionId = self.gaSession.sessionId
		self.dbSession.add(self.gaSession)		

	def saveMetrics(self, metrics):				
		self.gaSession.sessionMetrics = Metrics(self.gaSessionId, metrics)				
		
	def saveGenealogy(self, tree, history):
		self.gaSession.sessionGenealogy = Genealogy(self.gaSessionId, tree, history)

	def saveIndividual(self, gen, ind, fitness, genome, phenome):
		self.gaSession.sessionIndividuals.append(Individual(self.gaSessionId, gen, ind, fitness, genome, phenome))

	def whatSessions(self):
		gaSessions = self.dbSession.query(Sessions).all()	
		return [gaSession.sessionId for gaSession in gaSessions]	

	def loadSession(self, sessionId):
		gaSession = self.dbSession.query(Sessions).filter(Sessions.sessionId == sessionId).first()
		data = {}

		if gaSession.sessionMetrics:
			data['metrics'] = gaSession.sessionMetrics.metricsPickle

		if gaSession.sessionGenealogy:
			data['genealogy'] = {}
			data['genealogy']['tree'] = gaSession.sessionGenealogy.treePickle
			data['genealogy']['history'] = gaSession.sessionGenealogy.historyPickle

		if gaSession.sessionIndividuals:
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

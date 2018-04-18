from sqlalchemy import Table, create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, backref

import time

Base = declarative_base()

class Session(Base):
	__tablename__ ='sessions'

	pID = Column('id', Integer, primary_key=True)
	timestamp = Column('timestamp', String)	
	metrics = relationship('Metrics',uselist=False, back_populates='session')
	genealogy = relationship('Genealogy', uselist=False, back_populates='session')
	generations = relationship('Generation',back_populates='session')
	arguments = Column('arguments', String)

	def __init__(self, sessionTimeStamp,arguments):
		self.timestamp = sessionTimeStamp
		self.arguments = arguments

	def getIndividualsList(self):
		inds = []
		for gen in self.generations:
			for ind in gen.individuals:
				inds.append(ind)
		return inds

class Metrics(Base):
	__tablename__ = 'metrics'
	
	pID = Column('id', Integer, primary_key=True)
	session_id = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Session',uselist=False, back_populates='metrics')
	metricsPickle = Column('metrics_pickle', PickleType)

	def __init__(self, sessionId, metrics):
		self.metricsPickle = metrics		

class Genealogy(Base):
	__tablename__ = 'genealogy'

	pID = Column('id', Integer, primary_key=True)
	session_id = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Session',uselist=False, back_populates='genealogy')
	treePickle = Column('tree_pickle', PickleType)
	historyPickle = Column('history_pickle', PickleType)

	def __init__(self, sessionId, tree, history):
		self.treePickle = tree
		self.historyPickle = history


class Generation(Base):
	__tablename__ = 'generations'

	pID = Column('id', Integer, primary_key=True)
	genNumber = Column('gen_number',Integer)
	sessionId = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Session',uselist=False,back_populates='generations')
	novelGenes = relationship('Gene',back_populates='generation')
	individuals = relationship('Individual',back_populates='gen')

	def __init__(self, genNum):
		self.genNumber = genNum

association_table = Table('association_ind_gene', Base.metadata,
    Column('individual_id', Integer, ForeignKey('individuals.id')),
    Column('gene_id', Integer, ForeignKey('genes.id'))
)

class Gene(Base):
	__tablename__ = 'genes'

	pID = Column('id', Integer, primary_key=True)
	rawGene = Column('raw_gene',String)
	gen_id = Column(Integer, ForeignKey('generations.id'))
	generation = relationship('Generation',uselist=False)
	individuals = relationship('Individual',secondary=association_table,back_populates='genes')

	def __init__(self, rawGene):
		self.rawGene = str(rawGene).replace('[','').replace(']','').replace(',','').replace(' ','')


class Individual(Base):
	__tablename__ = 'individuals'

	pID = Column('id', Integer, primary_key=True)
	gen_id = Column(Integer, ForeignKey('generations.id'))
	gen = relationship('Generation',uselist=False)

	gh_index = Column('gh_index', Integer)

	fitness = Column('fitness', Numeric)
	genome = Column('genome',String)
	genomePickle = Column('genome_pickle', PickleType)
	phenomePickle = Column('phenome_pickle', PickleType)
	genes = relationship('Gene',secondary=association_table,back_populates='individuals')

	def __init__(self, individual, phenome):
		self.fitness = individual.fitness.values[-1]
		self.genome = str(individual).replace('[','').replace(']','').replace(',','').replace(' ','')
		self.genomePickle = individual
		self.phenomePickle = phenome
		self.gh_index = individual.history_index
		

	def addGene(self,gene):
		self.genes.append(gene)

	def getGenNumber(self):
		return self.gen.genNumber


class DatabaseConnection:

	def __init__(self, dbfile):		
		engine = create_engine('sqlite:///{}'.format(dbfile))						
		Base.metadata.create_all(bind=engine)
		dbSession = sessionmaker(bind=engine)
		self.dbSession = dbSession()
		
	def saveSession(self,arguments):
		self.gaSessionTimeStamp = time.strftime("%Y-%m-%d %H:%M:%S")
		self.gaSession = Session(self.gaSessionTimeStamp,arguments)
		self.gaSessionId = self.gaSession.pID
		self.dbSession.add(self.gaSession)		

	def saveMetrics(self, metrics):				
		self.gaSession.metrics = Metrics(self.gaSessionId, metrics)				
		
	def saveGenealogy(self, tree, history):
		self.gaSession.genealogy = Genealogy(self.gaSessionId, tree, history)

	def saveIndividual(self, ind):
		self.gaSession.sessionIndividuals.append(ind)

	def saveGeneration(self,gen):
		self.gaSession.generations.append(gen)

	def whatSessions(self):
		gaSessions = self.dbSession.query(Session).all()	
		return [gaSession.pID for gaSession in gaSessions]

	def getLastSession(self):
		ids = whatSessions()
		return getSession(ids[-1]) if len(ids) > 0 else None

	def getSession(self,sessionId):
		return self.dbSession.query(Session).filter(Session.pID == sessionId).first()

	#deprecated
	def loadSession(self, sessionId):
		gaSession = self.dbSession.query(Session).filter(Session.pID == sessionId).first()
		data = {}

		if gaSession.metrics:
			data['metrics'] = gaSession.metrics.metricsPickle

		if gaSession.genealogy:
			data['genealogy'] = {}
			data['genealogy']['tree'] = gaSession.genealogy.treePickle
			data['genealogy']['history'] = gaSession.genealogy.historyPickle

		if gaSession.generations:
			data['individuals'] = []
			for sessionGeneration in gaSession.generations:
				for sessionIndividual in sessionGeneration.individuals:
					individual = {}
			 		individual['gen'] = sessionGeneration.genNumber
			 		individual['gh'] = sessionIndividual.gh_index
			 		individual['fitness'] = sessionIndividual.fitness
					individual['genome'] = sessionIndividual.genomePickle
					individual['phenome'] = sessionIndividual.phenomePickle
					data['individuals'].append(individual)

		return data

	def commit(self):
		self.dbSession.commit()

	def close(self):
		self.dbSession.close()		

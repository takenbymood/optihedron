from sqlalchemy import Table, create_engine, Column, String, PickleType, Integer, Numeric, ForeignKey
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
	sessionGenerations = relationship('Generation',back_populates='session')
	arguments = Column('arguments', String)

	def __init__(self, sessionTimeStamp,arguments):
		self.sessionTimeStamp = sessionTimeStamp
		self.arguments = arguments

class Metrics(Base):
	__tablename__ = 'metrics'
	
	pID = Column('id', Integer, primary_key=True)
	session_id = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Sessions',uselist=False, back_populates='sessionMetrics')
	metricsPickle = Column('metrics_pickle', PickleType)

	def __init__(self, sessionId, metrics):
		self.metricsPickle = metrics		

class Genealogy(Base):
	__tablename__ = 'genealogy'

	pID = Column('id', Integer, primary_key=True)
	session_id = Column(Integer, ForeignKey('sessions.id'))
	session = relationship('Sessions',uselist=False, back_populates='sessionGenealogy')
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
	session = relationship('Sessions',uselist=False,back_populates='sessionGenerations')
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
	fitness = Column('fitness', Numeric)
	genome = Column('genome',String)
	genomePickle = Column('genome_pickle', PickleType)
	phenomePickle = Column('phenome_pickle', PickleType)
	genes = relationship('Gene',secondary=association_table,back_populates='individuals')

	def __init__(self, fitness, genome, phenome):
		self.fitness = fitness
		self.genome = str(genome).replace('[','').replace(']','').replace(',','').replace(' ','')
		self.genomePickle = genome
		self.phenomePickle = phenome
		

	def addGene(self,gene):
		self.genes.append(gene)


class DatabaseConnection:

	def __init__(self, dbfile):		
		engine = create_engine('sqlite:///{}'.format(dbfile))						
		Base.metadata.create_all(bind=engine)
		dbSession = sessionmaker(bind=engine)
		self.dbSession = dbSession()
		
	def saveSession(self,arguments):
		self.gaSessionTimeStamp = time.strftime("%Y-%m-%d %H:%M:%S")
		self.gaSession = Sessions(self.gaSessionTimeStamp,arguments)
		self.gaSessionId = self.gaSession.sessionId
		self.dbSession.add(self.gaSession)		

	def saveMetrics(self, metrics):				
		self.gaSession.sessionMetrics = Metrics(self.gaSessionId, metrics)				
		
	def saveGenealogy(self, tree, history):
		self.gaSession.sessionGenealogy = Genealogy(self.gaSessionId, tree, history)

	def saveIndividual(self, ind):
		self.gaSession.sessionIndividuals.append(ind)

	def saveGeneration(self,gen):
		self.gaSession.sessionGenerations.append(gen)

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

from sqlalchemy import create_engine, Column, String, Numeric, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship

import time

Base = declarative_base()

class Sessions(Base):
	__tablename__ ='sessions'

	sessionId = Column('sessionId', String, primary_key=True)	
	sessionMetrics = relationship('Metrics', backref='sessions')

	def __init__(self, sessionId):
		self.sessionId = sessionId

class Metrics(Base):
	__tablename__ = 'metrics'
	
	sessionId = Column(String, ForeignKey('sessions.sessionId'), primary_key=True)
	gen = Column('gen', Numeric, primary_key = True)
	fitnessStd = Column('std', Numeric) 
	fitnessMax = Column('max', Numeric)
	fitnessAvg = Column('avg', Numeric)
	fitnessMin = Column('min', Numeric)

	def __init__(self, sessionId, gen, fitnessStd, fitnessMax, fitnessAvg, fitnessMin):
		self.sessionId = sessionId
		self.gen = gen
		self.fitnessStd = fitnessStd
		self.fitnessMax = fitnessMax
		self.fitnessAvg = fitnessAvg
		self.fitnessMin = fitnessMin

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
		print metrics
		for metric in metrics:
			self.dbSession.add(Metrics(self.gaSessionId, metric['gen'],metric['std'],metric['max'],metric['avg'],metric['min']))			
		self.dbSession.commit()				
		

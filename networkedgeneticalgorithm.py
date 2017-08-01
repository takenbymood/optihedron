import array
import random
import math

import itertools
import operator

import numpy
from customMap import customMap
import networks

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import joblib

from joblib import Parallel, delayed, parallel_backend
from distributed.joblib import DaskDistributedBackend

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

def defaultAlgorithmEaSimple(pop,toolbox,stats,hof):
		return algorithms.eaSimple(pop,toolbox=toolbox,
		cxpb=0.5, mutpb=0.2, ngen=5,verbose=False,stats=stats,halloffame=hof)

def defaultTwoPoint(ind1, ind2):
    size = len(ind1)
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else: # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
        
    return ind1, ind2

def defaultEvalMax(individual):
    return sum(individual),

def defaultSelTournament(pop,k):
    return tools.selTournament(pop,k,3)

def defaultMutFlipBit(individual):
    return tools.mutFlipBit(individual,0.05)


class NetworkedGeneticAlgorithm:

    def __init__(
    	self,
    	genomeSize,
    	islePop,
    	evaluate = defaultEvalMax,
    	net = networks.createIslands(10),
    	mate = defaultTwoPoint,
    	mut = defaultMutFlipBit,
    	sel = defaultSelTournament,
    	subroutine = defaultAlgorithmEaSimple,
    	atMigration=lambda x: None):

        
        self.toolbox = base.Toolbox()

        # Attribute generator
        self.toolbox.register("attr_bool", random.randint, 0, 1)

        # Structure initializers
        self.toolbox.register("individual", tools.initRepeat, creator.Individual, self.toolbox.attr_bool, genomeSize)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

        self.toolbox.register("evaluate", evaluate)
        self.toolbox.register("mate", mate)
        self.toolbox.register("mutate", mut)
        self.toolbox.register("select", sel)
        self.toolbox.register("map", map)
        self.net = net
        self.subroutine = subroutine
        self.islePop = islePop
        self.atMigration = atMigration

    

    def genMetrics(self,gen,island,logbook):
        log = []
        for d in logbook:
            d['gen']+=gen
            d['island']=island
            log.append(d)
        return log

    def getTotalFitness(self,pop):
        return numpy.sum([p.fitness.values[-1] for p in pop])

    def accumulateStats(self,l):
        acc = []
        for key, group in itertools.groupby(l, lambda item: item["gen"]):
            glist = list(group)
            acc.append({'gen': key,
                'max':  numpy.max([item["max"] for item in glist])
                ,'min':  numpy.min([item["min"] for item in glist])
                ,'avg':    numpy.mean([item["avg"] for item in glist])
                ,'std':    math.sqrt(numpy.sum(map(lambda x: x*x,[item["std"] for item in glist])))
                })
        return acc


    def migration(self,islands):
        network = self.net
        rFitness = 1.0/numpy.sum(map(self.getTotalFitness,islands))
        s = random.uniform(0,1)
        t = 0.0
        for i in range(len(islands)):
            isle = islands[i]
            for ind in isle:
                f = ind.fitness.values[-1]*rFitness
                if f + t >= s:
                    edges = network.edges(network.nodes()[i])
                    if len(edges)>0:
                        edge = random.choice(edges)[-1]
                        dest = network.nodes().index(edge)
                        #print("migrating " + str(ind) + " from " + str(i)  + " to " + str(dest))
                        islands[dest][random.randint(0,len(islands[dest])-1)] = ind
                    return
                else:
                    t+=f

    def algorithm(self,pop):
        return self.subroutine(pop,self.toolbox,self.buildStats(),self.buildHOF())

    def buildStats(self):
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", numpy.mean)
        stats.register("std", numpy.std)
        stats.register("min", numpy.min)
        stats.register("max", numpy.max)
        return stats

    def buildHOF(self,size=1):
    	return tools.HallOfFame(size, similar=numpy.array_equal)
        

    def run(self,ngen,freq):
        self.metrics = []
        self.islands = [self.toolbox.population(n=self.islePop) for i in range(len(self.net))]
        for i in range(0, ngen, freq):
            results = self.toolbox.map(self.algorithm, self.islands)
            self.islands = [pop for pop, logbook in results]
            # print([d for d in logbook for pop, logbook in results])
            self.metrics += map(self.genMetrics,[i]*len(results),[n for n in range(len(results))], [logbook for pop, logbook in results])
            self.atMigration(self)
            self.migration(self.islands)
        self.metrics = [val for sublist in self.metrics for val in sublist]
        self.metrics = sorted(sorted(self.metrics, key=lambda k: k['island']), key=lambda k: k['gen']) 
        self.accMetrics = (self.accumulateStats(self.metrics))
        #self.metrics = list(self.accumulate(self.metrics))
        return self.islands, self.metrics, self.accMetrics






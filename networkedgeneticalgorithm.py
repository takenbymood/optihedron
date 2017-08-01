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


class NetworkedGeneticAlgorithm:

    def __init__(self,creator,genomeSize,evaluate,mate,mut,sel,net,subroutine,islePop,atMigration=lambda x: None):

        
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
        return (gen,island, numpy.max([d['max'] for d in logbook]))

    def getTotalFitness(self,pop):
        return numpy.sum([p.fitness.values[-1] for p in pop])

    def accumulateMetrics(self,lis):
        for k, g in itertools.groupby(lis, key=operator.itemgetter(0)):
            lis = [y for x in g for y in x[1:]]
            fits = [t[2] for t in lis]
            yield (k[0],numpy.mean(fits),numpy.std(fits),numpy.max(fits),numpy.min(fits))

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
    	return self.subroutine(pop,self.toolbox)
        

    def run(self,ngen,freq):
        self.metrics = []
        self.islands = [self.toolbox.population(n=self.islePop) for i in range(len(self.net))]
        for i in range(0, ngen, freq):
            results = self.toolbox.map(self.algorithm, self.islands)
            self.islands = [pop for pop, logbook in results]
            self.metrics.append(map(self.genMetrics,[i+freq]*len(results),[n for n in range(len(results))], [logbook for pop, logbook in results]))
            self.atMigration(self)
            self.migration(self.islands)
        self.metrics = list(self.accumulateMetrics(self.metrics))
        return self.islands, self.metrics






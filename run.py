import array
import random
import math

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


def cxTwoPointCopy(ind1, ind2):
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

def evalOneMax(individual):
    return sum(individual),

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

# Attribute generator
toolbox.register("attr_bool", random.randint, 0, 1)

# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, 100)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", cxTwoPointCopy)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
toolbox.register("map", map)

hof = tools.HallOfFame(1, similar=numpy.array_equal)
stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg", numpy.mean)
stats.register("std", numpy.std)
stats.register("min", numpy.min)
stats.register("max", numpy.max)

NSPOKES = 2
NISLES = 11
ISLESIZE = 5
CXPB=0.5
MUTPB=0.2
NGEN, FREQ = 100, 5
VERBOSE=False
mStar = networks.createMegaStar(NSPOKES,int(math.ceil((NISLES-3)*0.25)),int(math.floor((NISLES-3)*0.25)))

def algorithm(pop):
    return algorithms.eaSimple(pop,toolbox=toolbox,
        cxpb=CXPB, mutpb=MUTPB, ngen=FREQ,verbose=VERBOSE,stats=stats,halloffame=hof)

def saveMetrics(n,logbook):
    print(n, numpy.max([d['max'] for d in logbook]))

def getTotalFitness(pop):
    return numpy.sum([p.fitness.values[-1] for p in pop])

def migration(islands,network):
    rFitness = 1.0/numpy.sum(map(getTotalFitness,islands))
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
                    islands[dest][random.randint(0,len(islands[dest])-1)] = ind
                return
            else:
                t+=f


def main():
    random.seed(64)
    islands = [toolbox.population(n=ISLESIZE) for i in range(NISLES)]

    # Unregister unpicklable methods before sending the toolbox.
    toolbox.unregister("attr_bool")
    toolbox.unregister("individual")
    toolbox.unregister("population")

    for i in range(0, NGEN, FREQ):
        results = toolbox.map(algorithm, islands)
        islands = [pop for pop, logbook in results]
        map(saveMetrics,[n + i for n in range(len(results))], [logbook for pop, logbook in results])
        migration(islands,mStar)
    return islands

if __name__ == "__main__":
    main()
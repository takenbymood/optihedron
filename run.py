import array
import random
import math
import csv
import argparse

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

import networkedgeneticalgorithm as nga

NSPOKES = 2
NISLES = 11
ISLESIZE = 5
CXPB=0.5
MUTPB=0.2
NGEN, FREQ = 100, 2
VERBOSE=False
mStar = networks.createMegaStar(NSPOKES,int(math.ceil((NISLES-3)*0.25)),int(math.floor((NISLES-3)*0.25)))

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

def selTournament(pop,k):
    return tools.selTournament(pop,k,3)

def mutFlipBit(individual):
    return tools.mutFlipBit(individual,0.05)


def algorithm(pop,toolbox):
    hof = tools.HallOfFame(1, similar=numpy.array_equal)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)
    return algorithms.eaSimple(pop,toolbox=toolbox,
        cxpb=CXPB, mutpb=MUTPB, ngen=FREQ,verbose=VERBOSE,stats=stats,halloffame=hof)

def saveMetrics(lis):
    with open('ur file.csv','wb') as out:
        csv_out=csv.writer(out)
        csv_out.writerow(['gen','avg','std','min','max'])
        for row in lis:
            csv_out.writerow(row)


def main():

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

    random.seed(64)
    ga = nga.NetworkedGeneticAlgorithm(creator,
        100,
        evalOneMax,
        cxTwoPointCopy,
        mutFlipBit,
        selTournament,
        mStar,
        algorithm,
        10)
    results = ga.run(NGEN,FREQ)
    saveMetrics(results[-1])



if __name__ == "__main__":
    main()
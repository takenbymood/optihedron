import array
import random
import math
import csv
import argparse
import time

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

parser = argparse.ArgumentParser(description='')
parser.add_argument('-n','--ngen', default=100,type=int,
                    help='number of generations')
parser.add_argument('-f','--freq', default=2, type=int,
                    help='number of generations between migrations')
parser.add_argument('-d','--demes', default=10, type=int,
                    help='number of demes to run over')
parser.add_argument('-p','--pop', default=10, type=int,
                    help='population of each deme')
parser.add_argument('-gs','--genomesize', default=100, type=int,
                    help='number of bits in the genome')
parser.add_argument('-c','--cxpb', default=0.5,  type=float,
                    help='independant probability of crossover')
parser.add_argument('-m','--mutpb', default=0.2, type=float,
                    help='independant probability of mutation')
parser.add_argument('-mg','--migrations', default=1, type=int,
                    help='number of migrations to do each time')
parser.add_argument('-mpb','--mindpb', default=0.05, type=float,
                    help='independant probability of a bit flip')
parser.add_argument('-t','--tournsize', default=3, type=int,
                    help='size of selection tournaments')
parser.add_argument('-v','--verbose', default=False, action='store_true',
                    help='option to run in verbose mode')
parser.add_argument('-g','--graph', default='islands', 
                    choices=['singlet','islands','star','megastar'],
                    help='type of network to use')
parser.add_argument('-s','--seed', default=int(time.time()), type=int,
                    help='seed for the RNG')



args = parser.parse_args()

NSPOKES = 2
NISLES = args.demes
ISLESIZE = args.pop
CXPB=args.cxpb
MUTPB=args.mutpb
NGEN, FREQ = args.ngen, args.freq
VERBOSE=args.verbose
SEED = args.seed
TSIZE = args.tournsize
MINPDB = args.mindpb
GENOMESIZE = args.genomesize
MIGR = args.migrations

def saveMetrics(lis,filename='out.csv'):
    with open(filename,'wb') as out:
        csv_out=csv.DictWriter(out,lis[-1].keys())
        csv_out.writeheader()
        for row in lis:
            csv_out.writerow(row)

def evaluate(individual):
    return sum(individual),

def sel(pop,k):
    return tools.selTournament(pop,k,TSIZE)

def mut(individual):
    return tools.mutFlipBit(individual,MINPDB)

def algorithm(pop,toolbox,stats,hof):
    return algorithms.eaSimple(pop,toolbox=toolbox,
        cxpb=CXPB, mutpb=MUTPB, ngen=FREQ,verbose=VERBOSE,stats=stats,halloffame=hof)

def beforeMigration(ga):
    return

def afterMigration(ga):
    return

def main():

    if args.graph == 'singlet':
        network = networks.createSinglets(NISLES)
    elif args.graph == 'islands':
        network = networks.createIslands(NISLES)
    elif args.graph == 'star':
        network = networks.createStar(NISLES-1)
    elif args.graph == 'megastar':
        network = networks.createMegaStar(NSPOKES,int(math.ceil((NISLES-3)*0.25)),int(math.floor((NISLES-3)*0.25)))
    else:
        raw_input('malformed network option, continue with islands? (Enter)')

    random.seed(args.seed)
    ga = nga.NetworkedGeneticAlgorithm(
        genomeSize = GENOMESIZE,
        islePop = ISLESIZE,
        evaluate = evaluate,
        sel = sel,
        net = network,
        subroutine = algorithm,
        mut = mut,
        beforeMigration = beforeMigration,
        afterMigration = afterMigration)

    results = ga.run(NGEN,FREQ,MIGR)
    #print(results[0][0][0])
    saveMetrics(results[-1])



if __name__ == "__main__":
    main()
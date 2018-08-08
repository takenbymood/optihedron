import numpy
import random
from copy import deepcopy

#https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html

class IndArray(numpy.ndarray):

    def __init__(self, attributes):
        # Some initialisation with received values
        pass

    def __new__(subtype, shape, dtype=int, buffer=None, offset=0,
                strides=None, order=None, info=None):
        obj = super(IndArray, subtype).__new__(subtype, shape, dtype,
                                                buffer, offset, strides,
                                                order)
        obj.info = info
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', None)

    def __deepcopy__(self,memo):
        cls = self.__class__
        result = cls.__new__(cls,len(self))
        for i in range(len(self)):
            result[i] = deepcopy(self[i])
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        result.fitness = deepcopy(self.fitness)
        result.info = deepcopy(self.info)
        return result

def make_individual_valid(ind,params):
    ind.info = params
    return ind

def generate_individual(ind_class, size):
    individual = ind_class(size)
    for i in range(len(individual)):
        individual[i] = random.randint(0,1)
    individual = make_individual_valid(individual, 'test')
    return individual
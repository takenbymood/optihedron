import customPickle

from joblib import Parallel, delayed

def customMap(f, *iters):
    return Parallel(n_jobs=-1)(delayed(f)(*args) for args in zip(*iters))
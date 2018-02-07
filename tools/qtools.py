import subprocess

def runCmd(exe):
    p = subprocess.Popen(exe,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        retcode = p.poll()
        line = p.stdout.readline()
        yield line
        if retcode is not None:
            break

def hasRQJob(number):
    jobs = runCmd(['qstat',str(number).split('.')[0]])
    for line in jobs:
        columns = line.split()
        if columns[-2] in ('Q','R'): return True
    return False
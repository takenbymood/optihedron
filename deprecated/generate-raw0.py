from models import modelparticles as mpar
from tools import vectools
import parlammps
from nanoparticle import CoveredNanoParticlePhenome
from membranesimulation import MembraneSimulation
from tools import analysistools as atools
#from pathos import pools
from operator import itemgetter
import os
import pickle
import random
import numpy


patchy36 = mpar.modelpatchy36()
patchy42 = mpar.modelpatchy42()
patchy33 = mpar.modelpatchy33()
patchy32 = mpar.modelpatchy32()
patchy34 = mpar.modelpatchy34()
patchy30 = mpar.modelpatchy30()
liney36 = mpar.ptl(patchy36, 72)
liney42 = mpar.ptl(patchy42, 72)
liney33 = mpar.ptl(patchy33, 72)
liney32 = mpar.ptl(patchy32, 72)
liney34 = mpar.ptl(patchy34, 72)
liney30 = mpar.ptl(patchy30, 72)

epstotal = [250.0,240.0,230.0,220.0,210.0,200.0]

rots = []
rots.append([[0.0, 1.0, 0.0],0.0,0])
#rots.append([[0.0, 1.0, 0.0],3.141,0])
#rots.append([[-1.0,0.0, 0.0],1.571,0])
#rots.append([[1.0, 0.0, 0.0],1.571,0])
#rots.append([[0.0,-1.0, 0.0],1.571,0])
#rots.append([[0.0, 1.0, 0.0],1.571,0])


designmode = [patchy36, patchy42, patchy33, patchy32, patchy34, patchy30, liney36, liney42, liney33, liney32, liney34, liney30]
designname = ['patchy36', 'patchy42', 'patchy33', 'patchy32', 'patchy34', 'patchy30', 'liney36', 'liney42', 'liney33', 'liney32', 'liney34', 'liney30']


def evaluateNPWrapping(np,outFilename,runtime):    
    minFit = 1E-8
    noBud = False
    outHeaderSize = 9
    outData = {}

    nActiveLigands = 0
    npTotalEps = 0.0

    for l in np.ligands:
        nActiveLigands += 1
        npTotalEps += l.eps

    if(not os.path.exists(outFilename)):                                
            return minFit, noBud

    with open(outFilename, 'r+') as f:
        lines = f.readlines()
        ts = 0
        steps = []
        hPos = 0        
        for i in range(len(lines)):
            if str('ITEM: TIMESTEP') in lines[i]:
                hPos = i
                ts = int(lines[i+1])
                steps.append(ts)
                outData[ts]=[]                
            if(i-hPos>outHeaderSize):
                outData[ts].append(lines[i].replace("\n","").replace(" ",","))

    if len(outData[ts])<50:        
        return minFit, noBud 

    stepData = []
    budTime = False

    for s in steps:   
        outVectors = {}
        for line in outData[s]:
            slist = line.split(",")[1:]
            sId = line.split(",")[0]
            if(len(slist)<3):            
                return minFit, noBud
            if not int(slist[0]) in outVectors:
                outVectors[int(slist[0])] = []
            outVectors[int(slist[0])].append({'id':sId,'x':float(slist[1]),'y':float(slist[2]), 'z':float(slist[3]), 'c':int(slist[4])})

        cStep = []
        mStep = []
        boxsize = 20
        for key, value in outVectors.iteritems():
            budded = False
            for v in value:
                cIds = [c['id'] for c in cStep]
                if not v['c'] in cIds:
                    cStep.append({'id':v['c'],'size':1})
                else:
                    cId = 0
                    cCount = 0
                    for cI in cIds:
                        if cI == v['c']:
                            cId = cCount
                        cCount += 1
                    cStep[cId]['size'] += 1

            if key == 2:
                for v in value:
                    inrange = 0
                    fmag = 0
                    for v2 in outVectors[1]:
                        xd = v['x']-v2['x']
                        yd = v['y']-v2['y']
                        zd = v['z']-v2['z']
                        #squared magnitude of the difference
                        m = xd*xd+yd*yd+zd*zd                                          
                        if(m<25.0):
                            mStep.append(v2['id'])

            nLargeClusters = 0
            for v in sorted(cStep, key=itemgetter('size')):
                if v['size'] > 250:
                    nLargeClusters += 1
            budded = nLargeClusters > 1
                                       
            stepData.append({'timestep':s,'clusters':cStep,'magnitudes':mStep,'cNum':len(cStep),'mNum':len(mStep), 'budded': budded})

        if budded:
            if not budTime:
                budTime = s

    msum = stepData[-1]['mNum']

    if(msum == 0):        
        return minFit, noBud

    
    #reward = msum

    # penalty = PENALTYWEIGHT*(1.0-(float(npTotalEps)/(float(EPSMAX)*float(GENES))))*100 if float(EPSMAX)*float(nActiveLigands) > 0.0 else 0.0

    # reward = (float(BUDDINGREWARD) + float(penalty)) if stepData[-1]['budded'] else float(msum)
    BUDDINGREWARD = 400.0 
    reward = (float(BUDDINGREWARD)) if stepData[-1]['budded'] else float(msum)

    return reward, budTime


for indmode, indname in zip(designmode, designname):
    for rotdat in rots:
        for eps in epstotal:
            rotidx = rotdat[2]            
            RTIME = 25000
            TIMESTEP = 0.01    
            
            ind = indmode
            OUTDIR = 'models/raw/out{}'.format(rotidx)
            RUNDIR = 'models/raw/run{}'.format(rotidx)
            TEMPLATEDATAPATH = 'mem/template/data.template'
            TEMPLATEINPUTPATH = 'mem/template/in.template'
            EXPRPLACES = 1
            EPSPLACES = 0
            LIGEPS = float(eps)/float(numpy.sum(ind))
            LIGEPS = round(LIGEPS, 5)
            EPSMIN = LIGEPS
            EPSMAX = LIGEPS

            phenome = CoveredNanoParticlePhenome(ind,EXPRPLACES,EPSPLACES,EPSMIN,EPSMAX)            
            np = phenome.particle        
            sim = MembraneSimulation(
                '{}-{}-{}'.format(indname, eps, rotidx),
                np,
                RTIME,
                TIMESTEP,            
                OUTDIR,
                RUNDIR,            
                TEMPLATEDATAPATH,
                TEMPLATEINPUTPATH,
                rAxis = rotdat[0],
                rAmount = rotdat[1]
                )
            ScriptPath = os.path.join(sim.filedir,sim.scriptName)
            sim.saveFiles()    
            parlammps.runSimSerial(ScriptPath)
            outFilePath = os.path.join(sim.outdir,sim.outName)
            sim.postProcessOutput(outFilePath)
            outFilePathProcessed = outFilePath+'a'
            f = 1E-8
            b = False
            f,b = evaluateNPWrapping(np,outFilePath,RTIME)
            sim.postProcessOutput(outFilePath)
            contactData, maxcontactligs = atools.measureLigandContact(outFilePathProcessed)
            dumpdata = [contactData, maxcontactligs, np, rotdat[0], rotdat[1], f, b]
            pickleName = '{}-{}-{}'.format(indname, eps, rotidx)            
            pickle.dump(dumpdata, open('models/raw/pickles{}/{}.p'.format(rotidx,pickleName), 'wb'))
            if os.path.exists(outFilePathProcessed):
                os.remove(outFilePathProcessed)
            if os.path.exists(outFilePath):
                os.remove(outFilePath)
            sim.deleteFiles()                            
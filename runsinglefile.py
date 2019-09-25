import membranesimulation as mb
import argparse
import os
import pickle
import nanoparticle
import parlammps


parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', default='', type=str, 
                    help='input lammps script file')
parser.add_argument('-o','--out', default='', type=str, 
                    help='output directory')

args = parser.parse_args()

istr = "010000011000100010101000101000001100000010010111100001101001110000010000"
individual = [int(i) for i in istr]

particle = nanoparticle.CoveredNanoParticlePhenome(individual,1,0,10,10)

np = particle.particle

sim = mb.MembraneSimulation(
        args.input.split('/')[-1].split('.')[0],
        np,
        25000,
        0.01,        
        args.out,
        os.path.dirname(args.input),
        "/Users/joelforster/Projects/optihedron/mem/template/data.template",
        "/Users/joelforster/Projects/optihedron/mem/template/in.template",
        rAxis=(0,1,0),
        rAmount=0        
        )

sim.saveFiles()
scriptPath = os.path.join(sim.filedir,sim.scriptName)
parlammps.runSimSerial(scriptPath)
outFilePath = os.path.join(sim.outdir,sim.outName)
sim.postProcessOutput(outFilePath)

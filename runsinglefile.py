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

istr = "111111111111111111110011100110011011111011111111001111101110011111101101"
individual = [int(i) for i in istr]

particle = nanoparticle.CoveredNanoParticlePhenome(individual,1,0,10,10)

np = particle.particle

sim = mb.MembraneSimulation(
        args.input.split('/')[-1].split('.')[0],
        np,
        5000,
        0.01,        
        args.out,
        os.path.dirname(args.input),
        "/Users/joelforster/Projects/optihedron/mem/template/data.template",
        "/Users/joelforster/Projects/optihedron/mem/template/in.template",
        rAxis=(1,1,0),
        rAmount=3.141        
        )

sim.saveFiles()
scriptPath = os.path.join(sim.filedir,sim.scriptName)
parlammps.runSimSerial(scriptPath)
outFilePath = os.path.join(sim.outdir,sim.outName)
sim.postProcessOutput(outFilePath)

import ovito
import argparse
import sys
import math
import os

from ovito.vis import Viewport, RenderSettings
from ovito import dataset
from ovito.modifiers import *

parser = argparse.ArgumentParser()

parser.add_argument('-t','--test', default=False, action='store_true',
                    help='option to split into directories based on some test')

argPathGroup = parser.add_mutually_exclusive_group()
argPathGroup.add_argument('-dir','--directory', default='', type=str, 
                    help='path to the directory containing the files to render')
argPathGroup.add_argument('-f','--file', default='', type=str, 
                    help='path to the file to render')

parser.add_argument('-o','--out', default=os.path.dirname(os.path.realpath(__file__)), type=str, 
                    help='path to the output directory')

args = parser.parse_args()

files = []

if args.file != '' and '.xyza' in args.file:
    # single file render
    files.append(args.file)
elif args.directory != '':
    # multiple files
    for f in os.listdir(args.directory):
        if '.xyza' in f:
            files.append(os.path.join(args.directory,f))

def nanoparticleModifiers(node):
    node.modifiers.append(SelectExpressionModifier(expression = '(affinity == 0) && ParticleType > 2'))
    radMod = ComputePropertyModifier()
    radMod.output_property = "Radius"
    radMod.expressions = ["Selection ? 0.0001 : ParticleType >2 ? 0.75 : ParticleType != 1 ? 4.0 : 0.5"]
    node.modifiers.append(radMod)
    modifier = ColorCodingModifier(
    particle_property = "affinity",
    gradient = ColorCodingModifier.Hot(),
    start_value=15,
    end_value=0
    )
    node.modifiers.append(modifier)
    node.compute()
    return


def renderScene(inpath,outpath,vp,timestep=0,size=(800,600),modifiers=[],endframe=False):
    # file imports
    node = ovito.io.import_file(inpath, multiple_frames = True)
    node.add_to_scene()
    node.source.cell.display.render_cell = False
    # viewport setup
    vp.zoom_all()
    # animation properties
    if not endframe:
        dataset.anim.current_frame=timestep
    else:
        dataset.anim.current_frame=dataset.anim.last_frame
    # modifier function injection
    for mod in modifiers:
        mod(node)
    # rendering
    vp.render(RenderSettings(filename=outpath, size=size,generate_alpha=True))
    node.remove_from_scene()
    return vp

for f in files:
    simName = os.path.basename(os.path.normpath(f)).split('.')[0].split('_')[0]

    
    # nanoparticle specific stuff
    if args.test:
        simName = os.path.join('nonbudding',simName) if '-1' in f else os.path.join('budding',simName)
    # end of nanoparticle specific stuff!


    # outdir = os.path.join(args.out,simName) if args.test else args.out
    outdir = args.out
    filePath = os.path.join(outdir,simName+".png")

    renderScene(f,filePath,Viewport(type = Viewport.Type.TOP),timestep=100,modifiers=[nanoparticleModifiers],endframe=True)
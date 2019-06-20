import ovito
import argparse
import sys
import math
import os
import time

from ovito.vis import Viewport, RenderSettings
from ovito import dataset
from ovito.modifiers import *
from matplotlib import pyplot as plt
import gc 


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

parser.add_argument('-l','--limit', default=-1, type=int, 
                    help='maximum number of scenes to render')

parser.add_argument('-g','--grid', default=False, action='store_true',
                    help='option to produce a 3x3 grid of images')

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

    transparencyModifier = ComputePropertyModifier()
    transparencyModifier.output_property = "Transparency"
    transparencyModifier.expressions = ["ParticleType == 1 ? Position.Z < -3 ? 0.85 : 0.8 : 0"]
    node.modifiers.append(transparencyModifier)

    node.compute()
    return

single = not args.grid 

def renderScene(inpath,outpath,timestep=0,size=(1024,720),modifiers=[],endframe=False):
    # file imports
    node = ovito.io.import_file(inpath, multiple_frames = True)
    node.add_to_scene()

    #ortho viewport setup
    vp = Viewport(type = Viewport.Type.ORTHO)
    vp.fov=30
    vp.camera_pos = (5, -5, 5)
    vp.camera_dir = (-1, 1, -0.8)
    node.source.cell.display.render_cell = False

    # top viewport setup
    #vp = Viewport(type = Viewport.Type.TOP)
    #vp.zoom_all()

    # animation properties
    if not endframe:
        dataset.anim.current_frame=timestep
    else:
        dataset.anim.current_frame=dataset.anim.last_frame
    # modifier function injection
    for mod in modifiers:
        mod(node)
    # rendering
    img = vp.render(RenderSettings(filename=outpath,size=size,generate_alpha=True))
    node.remove_from_scene()
    del node, vp, modifiers,img
    gc.collect()
    return

counter = 0

outdir = args.out


for f in files:

    simName = os.path.basename(os.path.normpath(f)).split('.')[0].replace('baseball','bb')

    
    # nanoparticle specific stuff
    if args.test:
        simName = os.path.join('nonbudding',simName) if '-1' in f else os.path.join('budding',simName)
    # end of nanoparticle specific stuff!



    # outdir = os.path.join(args.out,simName) if args.test else args.out
    

    if os.path.exists(os.path.join(outdir,simName+"_grid.png")):
        continue
    if counter >= args.limit:
        print('limit reached, exiting')
        break
    if single:
        filePath = os.path.join(outdir,simName+"_ortho.png")
        renderScene(f,filePath,timestep=0,modifiers=[nanoparticleModifiers],endframe=True)
    else:
        interval = 31.25
        filePaths=[]
        for i in range(9):
            ts = int(i*interval)
            if ts > 250:
                ts = 250
            filePath = os.path.join(outdir,simName+"_ortho"+"_"+str(i)+".png")
            filePaths.append(filePath)
            renderScene(f,filePath,timestep=ts,modifiers=[nanoparticleModifiers],endframe=False)
        plt.axis('off')
        plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
        fig = plt.figure(figsize = (16,9))
        
        imnum = 0
        for i,p in enumerate(filePaths):
            img = plt.imread(p)
            ax = fig.add_subplot(3,3,i+1)
            ax.axis('off')
            ax.set_title(int(i*interval)*100)
            ax.imshow(img)
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(os.path.join(outdir,simName+"_grid.png"))
        fig.clf()
        for i,p in enumerate(filePaths):
            if os.path.exists(p):
                os.remove(p)
        counter += 1
        plt.close()
        gc.collect()
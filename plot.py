import numpy as np
import matplotlib.pyplot as plt
import math
import os
import argparse

parser = argparse.ArgumentParser(description='Plot fitness logs from GA')
parser.add_argument('--in','-i', dest='filepath', required=False,
                  	default=None,
                    help='csv file for plotting')

parser.add_argument('--add','-a', action='append', dest='addfilepath', 
                    required=False, default=[],
                    help='additional csv files for plotting')
                    
parser.add_argument('--out','-o', dest='output', required=False,
                  	default=None,
                    help='output filename')

parser.add_argument('--figx', '-x', dest='figurex', required=False,
                     default=10, type=int,
                     help='output figure x dimension')

parser.add_argument('--figy', '-y', dest='figurey', required=False,
                     default=10, type=int,
                     help='output figure y dimension')
                    
plt.rcParams.update({'font.size': 14})

args = parser.parse_args()

outfile = args.output if args.output != None and args.output != "" else "plot.png"

def tworound(x, base=2):
    return int(base * round(float(x)/base))

def datafromtxt(filepath):
    return np.genfromtxt(filepath, delimiter=',', skip_header=1,
                     skip_footer=0, names=['STD','MAX','AVG','GEN','MIN'])

# example data
data = []
fname = []
data.append(datafromtxt(args.filepath))
fname.append(os.path.basename(args.filepath))
for afilepath in args.addfilepath:
    data.append(datafromtxt(afilepath))
    fname.append(os.path.basename(afilepath))

                     # These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)    

x = []
yorig = []
y = []
y2 = []
err = []
#errN = map(math.sqrt,data['NV'])
stdErr = []
errPerc = []
yerr = []
for datum in data:
    x.append(datum['GEN'])
    yorig.append(datum['AVG'])
    y.append(datum['AVG'])
    y2.append(datum['MAX'])
    err.append(datum['STD'])
    stdErr.append(datum['STD']/10)
    errPerc.append((datum['STD']/10)/datum['AVG'])
    yerr.append(datum['AVG']*((datum['STD']/10)/datum['AVG']))

fig = plt.figure(figsize=(args.figurex,args.figurey))

ax1 = fig.add_subplot(111)
ax1.spines["top"].set_visible(False)  
ax1.spines["right"].set_visible(False)  

#ax1.set_title("Genome Fitness")    
ax1.set_xlabel('Generation Number')
ax1.set_ylabel('Fitness')

#plt.fill_between(x, y-yerr,y+yerr, color="#3F5D7D")  
for n in range(0,100,2):
    plt.plot([-1,np.max(x)+1], [n,n], "--", lw=0.5, color="black", alpha=0.3)  
    
plt.tick_params(axis="both", which="both", bottom="off", top="off",labelbottom="on", left="off", right="off", labelleft="on") 

for x_i, y_i, y2_i, yerr_i, fname_i, color_i in zip(x, y, y2, yerr, fname, tableau20[::2]):
    ax1.plot(x_i,y2_i,lw=1.5,label='Maximum {}'.format(fname_i),color=color_i)
    ax1.errorbar(x_i,y_i, yerr=yerr_i, markersize='3.5',capsize=2.5, fmt='o-', label='Average {}'.format(fname_i), color=color_i)    

plt.legend()


plt.ylim(tworound(np.min(y)-3),tworound(np.max(y2)+2))

#plt.ylim(20,60)
plt.xlim(0,np.max(x)+1)

plt.draw()
plt.savefig(outfile)
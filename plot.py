import numpy as np
import matplotlib.pyplot as plt
import math
import argparse

parser = argparse.ArgumentParser(description='Plot fitness logs from GA')
parser.add_argument('--in','-i', dest='filepath', required=False,
                  	default=None,
                    help='csv file for plotting')
                    
parser.add_argument('--out','-o', dest='output', required=False,
                  	default=None,
                    help='output filename')
                    
plt.rcParams.update({'font.size': 14})

args = parser.parse_args()

outfile = args.output if args.output != None and args.output != "" else "plot.png"

def tworound(x, base=2):
    return int(base * round(float(x)/base))

# example data
data = np.genfromtxt(args.filepath, delimiter=',', skip_header=1,
                     skip_footer=0, names=['STD','MAX','AVG','GEN','MIN'])

                     
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

x = data['GEN']
yorig = data['AVG']
y = yorig
y2 = data['MAX']
err = data['STD']
#errN = map(math.sqrt,data['NV'])
stdErr = err/10
errPerc = (stdErr/yorig)
yerr = y*errPerc

fig = plt.figure()

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

ax1.plot(x,y2,color=tableau20[2],lw=1.5,label='Maximum')
ax1.errorbar(x,y, yerr=yerr,color=tableau20[0], markersize='3.5',capsize=2.5, fmt='o')
ax1.plot(x, y, color=tableau20[0],lw=2, label='Average')

plt.legend(bbox_to_anchor=(0.675, 0.2), loc=2, borderaxespad=0.)


plt.ylim(tworound(np.min(y)-3),tworound(np.max(y2)+2))

#plt.ylim(20,60)
plt.xlim(0,np.max(x)+1)


plt.draw()
plt.savefig(outfile)
#plt.show()
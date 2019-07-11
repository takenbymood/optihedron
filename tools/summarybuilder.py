import analysistools as atools
import pandas as pd
import pickle
import os

import argparse

parser = argparse.ArgumentParser(description='')

parser.add_argument('-o','--output', default="", type=str, 
                    help='output directory for the generated files')
parser.add_argument('-i','--input', default='', type=str, 
                    help='input directory, files must be in xyza format')

args = parser.parse_args()

xyzaPath = args.input
outPath = args.output

s = {}

try:
    s = atools.generateSummaries(xyzaPath)
except:
    print('something went wrong')

with open(os.path.join(outPath,'trajectories.pickle'), 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(s, f, pickle.HIGHEST_PROTOCOL)

df = pd.DataFrame()
data = []
bmsd = []
nmsd = []
for k,v in s.iteritems():
    if (len(bmsd) > 0 and len(nmsd) > 0) and (len(v['msd']) != len(bmsd) or len(v['msd']) != len(nmsd)):
        continue
    data.append((k,v['density'],v['clustering'],v['bt'],'budding' if v['bt'] > 0.0 else 'non budding'))
    try:
        if v['bt'] > 0.0:
            bmsd = [bmsd[i]+t for i,t in enumerate(v['msd'])] if len(bmsd) > 0 else [t for i,t in enumerate(v['msd'])]
        else:
            nmsd = [nmsd[i]+t for i,t in enumerate(v['msd'])] if len(nmsd) > 0 else [t for i,t in enumerate(v['msd'])]
    except:
        print('something went wrong with msd')

try:
    if len(bmsd) > 0:
        bmsd = [t/float(len(bmsd)) for t in bmsd]
    if len(nmsd) > 0:
        nmsd = [t/float(len(nmsd)) for t in nmsd]
except:
    print('couldnt get msd avg')

print bmsd
print nmsd

df = pd.DataFrame(data, columns = ['file','den','cls','bud', 'cat' ]) 
df.to_csv(os.path.join(outPath,'trajectories-summary.csv'),index=False)
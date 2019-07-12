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
for k,v in s.iteritems():
    data.append((k,v['density'],v['clustering'],v['bt'],'budding' if v['bt'] > 0.0 else 'non budding'))

df = pd.DataFrame(data, columns = ['file','den','cls','bud', 'cat' ]) 
df.to_csv(os.path.join(outPath,'trajectories-summary.csv'),index=False)
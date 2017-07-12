#!/usr/bin/env python

# parse nexus output such as:
#   hse-0/dmc/hse-0-dmc  series 0  -10.504233 +/- 0.003069    1.0   0.481375 +/- 0.022337    1.0   0.0458 

def exx_from_myid(myid):
  exx = int( myid.split('/')[0].replace('hse-','') )
  return exx

if __name__ == '__main__':
  import os
  import pandas as pd

  col_map = {'energy_mean':3,'energy_error':5,'variance_mean':5,'variance_error':7,'iqmc':2}
  data = {'exx':[],'energy_mean':[],'energy_error':[],'variance_mean':[],'variance_error':[],'iqmc':[]}
  
  json_fname = 'qmc.json'
  if not os.path.isfile(json_fname):
    with open('output.txt','r') as f:
      for line in f:
        tokens = line.split()
        if len(tokens) == 12:
          myid = tokens[0]
          exx = exx_from_myid(myid)

          data['exx'].append(exx)
          for obs in col_map.keys():
            data[obs].append(tokens[col_map[obs]])
          # end for
        # end if
      # end for
    # end with

    df = pd.DataFrame(data)
    df.to_json(json_fname)
  else :
    df = pd.read_json(json_fname)
  # end if

  df.sort_values('exx',inplace=True)
  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  for iqmc in df['iqmc'].unique():
    mydf = df[df['iqmc']==iqmc].sort_values('exx')
    ax.errorbar(mydf['exx'],mydf['energy_mean'],yerr=mydf['energy_error'].values)
  plt.show()

#!/usr/bin/env python

def collect_dft(json_fname):
  import subprocess as sp
  proc = sp.Popen(['find','../','-path','*-scf.out'],stdout=sp.PIPE,stderr=sp.PIPE)
  out,err = proc.communicate()

  import qe_reader as qer
  data = {'volume':[],'energy':[],'pseudo':[]}
  for fname in out.split('\n')[:-1]:
    # fname e.g.: ../bfd/pbe-102/scf/pbe-102-scf.out
    tokens = fname.split('/')
    pseudo = tokens[-4]
    #func,strain = tokens[3]; float(strain)/100.
    
    energy = qer.read_first_energy(fname)
    from mmap import mmap
    with open(fname,'r+') as f:
      mm = mmap(f.fileno(),0)
    # end with
    idx = mm.find('unit-cell volume')
    mm.seek(idx)
    line = mm.readline()
    vol = float(line.split('=')[1].split()[0])

    data['volume'].append(vol)
    data['energy'].append(energy)
    data['pseudo'].append(pseudo)
  # end for

  import pandas as pd
  df = pd.DataFrame(data)
  df.to_json(json_fname)
  return df
# end def collect_dft

if __name__ == '__main__':
  import os
  import pandas as pd
  fname = 'cold_curve.json'
  if not os.path.isfile(fname):
    df = collect_dft(fname)
  else:
    df = pd.read_json(fname)
  # end if

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel(r'volume (bohr$^3$)',fontsize=16)
  ax.set_ylabel('total energy - min (Ry)',fontsize=16)
  for pseudo in df['pseudo'].unique():
    mydf = df[ df['pseudo']==pseudo ].sort_values('volume')
    ax.plot(mydf['volume'],mydf['energy']-mydf['energy'].min(),'o-',label=pseudo)
  # end for
  plt.show()
# end __main__

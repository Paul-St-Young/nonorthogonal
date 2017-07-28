#!/usr/bin/env python
import pandas as pd

def ldau_from_myid(myid):
  # myid e.g. u0.50/dmc/lda-0-ecut320-dmc
  tokens = myid.split('/')
  subdir = tokens[-3]
  myu = float(subdir.replace('u',''))
  return pd.Series({'ldau':myu})
# end def

if __name__ == '__main__':
  import sys
  sys.path.insert(0,'../../../utils')
  from parsing import parse_qmcas_output

  result_fname = '../scan/ldau/lda-0-ecut320/results.dat'
  df0 = parse_qmcas_output(result_fname)
  settings = df0['myid'].apply(ldau_from_myid)
  df = pd.concat([settings,df0],axis=1)

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlim(-0.5,10.5)
  ax.set_xlabel('Hubbard U (eV)',fontsize=16)
  ax.set_ylabel('Total energy (Ha)',fontsize=16)

  for method in df.method.unique():
    sel = (df.method==method)
    if method == 'DMC':
      sel = (df.method==method) & (df.iqmc==2)
    # end if
    mydf = df.loc[sel].sort_values('ldau')

    ax.errorbar(mydf['ldau'],mydf['energy_mean'],yerr=mydf['energy_error'],label=method)
  # end for

  ax.legend(loc='upper left')
  fig.savefig('ldau_scan.png',dpi=320)
  plt.show()


# end __main__

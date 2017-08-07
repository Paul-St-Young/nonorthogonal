#!/usr/bin/env python
import numpy as np
import pandas as pd

def exx_from_myid(myid):
  # myid e.g. hse-50-ecut320/dmc/hse-50-ecut320-dmc
  tokens = myid.split('/')
  subdir = tokens[-3]
  exx = float(subdir.split('-')[1])
  return pd.Series({'exx':exx})
# end def

def best_improvement(myx,myy,myye):
  min_idx = np.argmin(myy)
  max_idx = np.argmax(myy)

  diff_mean = myy[max_idx] - myy[min_idx]
  diff_error= np.sqrt( myye[max_idx]**2. + myye[min_idx]**2. )
  return min_idx,myx[min_idx],diff_mean,diff_error
# end def

if __name__ == '__main__':
  nMn   = 2 # !!!! hard code the number of Mn atoms
  ha2ev = 27.21138602
  eunit = ha2ev/nMn # eV/f.u.
  #eunit = 1. # ha

  import sys
  sys.path.insert(0,'../../../utils')
  from parsing import parse_qmcas_output

  result_fname = '../cont_scan/exx/results.dat'
  df0 = parse_qmcas_output(result_fname)
  settings = df0['myid'].apply(exx_from_myid)
  df = pd.concat([settings,df0],axis=1)

  import matplotlib.pyplot as plt
  from matplotlib.ticker import FormatStrFormatter

  fig,ax = plt.subplots(1,1)
  ax.set_xlim(-1,101)
  ax.set_xlabel('HSE exact exachange fraction',fontsize=16)
  ax.set_ylabel('Total energy/f.u. (eV)',fontsize=16)

  for method in df.method.unique():
    sel = (df.method==method)
    if method == 'DMC':
      sel = (df.method==method) & (df.iqmc==2)
    # end if
    mydf = df.loc[sel].sort_values('exx')
    myx  = mydf['exx'].values
    myy  = mydf['energy_mean'].values*eunit
    myye = mydf['energy_error'].values*eunit

    idx,xmin,ydiff_mean,ydiff_error = best_improvement(myx,myy,myye)
    print xmin, ydiff_mean, ydiff_error

    ax.errorbar(myx,myy,yerr=myye,label=method)
  # end for

  ax.get_yaxis().set_major_formatter( FormatStrFormatter("%7.2f") )
  ax.legend(loc='upper left')
  fig.tight_layout()
  fig.savefig('qmc_exx_scan.png',dpi=320)
  plt.show()

# end __main__

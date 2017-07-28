#!/usr/bin/env python
import pandas as pd

def exx_from_myid(myid):
  # myid e.g. hse-50-ecut320/dmc/hse-50-ecut320-dmc
  tokens = myid.split('/')
  subdir = tokens[-3]
  exx = float(subdir.split('-')[1])
  return pd.Series({'exx':exx})
# end def

if __name__ == '__main__':
  import sys
  sys.path.insert(0,'../../../utils')
  from parsing import parse_qmcas_output

  df0 = parse_qmcas_output('../scan/exx/results.dat')
  settings = df0['myid'].apply(exx_from_myid)
  df = pd.concat([settings,df0],axis=1)


  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlim(-1,101)
  ax.set_xlabel('HSE exact exachange fraction',fontsize=16)
  ax.set_ylabel('Total energy (Ha)',fontsize=16)

  for method in df.method.unique():
    sel = (df.method==method)
    if method == 'DMC':
      sel = (df.method==method) & (df.iqmc==2)
    # end if
    mydf = df.loc[sel].sort_values('exx')

    ax.errorbar(mydf['exx'],mydf['energy_mean'],yerr=mydf['energy_error'],label=method)
  # end for

  ax.legend()
  plt.show()


# end __main__

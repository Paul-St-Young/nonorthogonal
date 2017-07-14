#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_styles():
  color_map_basis = {
  'double-zeta':'k',
  'triple-zeta':'g',
  'quadruple-zeta':'b'
  }

  ls_map_wf = {
    'Slater':':',
    'Slater-Jastrow':'-'
  }

  marker_map_method = {
    'phfmol':'o',
    'VMC':'^',
    'DMC':'s'
  }

  return color_map_basis, ls_map_wf, marker_map_method
# end def plot_styles

if __name__ == '__main__':

  color_map_basis, ls_map_wf, marker_map_method = plot_styles()

  df = pd.read_json('energy.json')

  line_map = {
    'line':[],
    'method':[],
    'wf':[],
    'basis':[]
  }

  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('number of determinants',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)
  ax.set_xlim(0,51)

  for method in df.method.unique():
    for wf in df.wf.unique():
      for basis in df.basis.unique():
        if method == 'DMC':
          sel = (df['method']==method) & (df['wf']==wf) & (df['basis']==basis)\
                  & (df['iqmc']==2)
        else:
          sel = (df['method']==method) & (df['wf']==wf) & (df['basis']==basis)
        # end if

        # get data
        mydf = df.loc[sel].sort_values('ndet')
        myx  = mydf['ndet'].values
        myy  = mydf['energy_mean'].values
        myye = mydf['energy_error'].values

        # get plotting style
        marker_size = 5
        if method != 'phfmol':
          marker_size *= 2
        mycolor = color_map_basis[basis]

        # plot
        line = ax.plot(myx,myy
          ,c      = mycolor
          ,ls     = ls_map_wf[wf]
          ,marker = marker_map_method[method]
          ,ms     = marker_size
          ,fillstyle = 'none'
        )
        ax.errorbar(myx,myy,yerr=myye,fmt='none',ecolor=mycolor)

        # save reference
        for name,obj in zip(['line','method','wf','basis'],[line[0],method,wf,basis]):
          line_map[name].append(obj)
        # end for

      # end for basis
    # end for wf
  # end for method

  # make legend 
  loc_map = {
   'method':'upper center'
  ,'wf':(0.65,0.65)
  ,'basis':'upper right'
  }
  for name in ['method','wf','basis']:
    mylines = []
    labels  = []
    for key in df[name].unique():
      idx = line_map[name].index(key)
      labels.append(key)
      mylines.append(line_map['line'][idx])
    # end for method
    leg = ax.legend(handles=mylines,labels=labels,loc=loc_map[name])
    ax.add_artist(leg)
  # end for

  fig.tight_layout()
  plt.show()

# end if

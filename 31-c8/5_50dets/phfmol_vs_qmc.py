#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def ndet_from_myid(myid):
  # e.g. myid = 'mdets/detsci24_af/c2'
  tokens = myid.split('/')
  det_str = tokens[-2]
  ndet_str = det_str.replace('detsci','').replace('_af','')
  return int(ndet_str)
# end def

if __name__ == '__main__':
  import sys
  sys.path.insert(0,'../../utils')
  from parsing import parse_qmcas_output

  # re-format phfmol data
  import pandas as pd
  phfdat = np.loadtxt('gen_dets/ewald.dat')
  pdf = pd.DataFrame(phfdat,columns=['ndet','energy_mean'])
  pdf['energy_error'] = 0.0
  pdf['method'] = 'phfmol'
  pdf['wf'] = 'Slater'

  # read QMC data
  df0 = parse_qmcas_output('s_vmc.dat')
  df0['wf'] = 'Slater'
  df0['ndet'] = df0['myid'].apply(ndet_from_myid)
  df1 = parse_qmcas_output('sj_qmc.dat')
  df1['wf'] = 'Slater-Jastrow'
  df1['ndet'] = df1['myid'].apply(ndet_from_myid)

  # combine all dataframes
  df = pd.concat([pdf,df0,df1]).reset_index(drop=True)
  df['basis'] = 'double-zeta'
  df['ndet'] += 1

  import sys
  sys.path.insert(0,'../../29-c2-basis')
  from show_energy import plot_styles
  color_map_basis, ls_map_wf, marker_map_method = plot_styles()

  # plot
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('number of determinants',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)
  ax.set_xlim(0,51)

  line_map = {
    'line':[],
    'method':[],
    'wf':[],
    'basis':[]
  }

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
        # ed for
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
  fig.savefig('c2-8-dz.png',dpi=300)
  plt.show()

#!/usr/bin/env python

def exx_from_myid(myid):
  # e.g. myid = 'hse-11/dmc/hse-11-dmc'
  tokens = myid.split('/')
  subdir = tokens[-3]
  exx_str = subdir.split('-')[1]
  return float(exx_str)

if __name__ == '__main__':
  import sys
  sys.path.insert(0,'../../utils')
  from parsing import parse_qmcas_output
  from plotting import plot_styles
  color_map_basis, ls_map_wf, marker_map_method = plot_styles()

  df = parse_qmcas_output('default/dmc.dat')
  df['exx'] = df['myid'].apply(exx_from_myid)
  df.sort_values('exx',inplace=True)

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('HSE exact exchange fraction (%)',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)

  for method in df.method.unique():
    if method == 'DMC':
      sel = (df.method == method) & (df.iqmc == 2)
    else:
      sel = (df.method == method)
    # end if

    myx = df.loc[sel,'exx'].values
    myy = df.loc[sel,'energy_mean'].values
    myye= df.loc[sel,'energy_error'].values
      
    ax.errorbar(myx,myy,yerr=myye,label=method,c='r' # PW
      ,marker=marker_map_method[method]
      ,ls=ls_map_wf['Slater-Jastrow']# all SJ wf
      ,fillstyle='none'
      ,markersize=10
    ) 
  # end for

  ax.set_ylim(-10.601,-10.25)
  ax.legend()
  fig.tight_layout()
  fig.savefig('c2-exx.png',dpi=480)
  plt.show()
# end __main__

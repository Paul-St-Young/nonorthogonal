#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def func_ecut_exx_from_input(inp_dict):
  func = inp_dict['system']['input_dft']
  ecut = float(inp_dict['system']['ecutwfc'])
  exx  = float(inp_dict['system']['exx_fraction'])
  return pd.Series({'func':func,'ecut':ecut,'exx':exx})
# end def

if __name__ == '__main__':
  df0 = pd.read_json('../mno4_qe_ecut_exx.json')

  # get func, ecut and exx
  settings = df0['input'].apply(func_ecut_exx_from_input)
  df = pd.concat([settings,df0],axis=1)

  fig   = plt.figure()

  label_fs = 14
  # energy vs. ecut
  iplot = 1
  ax = fig.add_subplot(2,2,iplot)
  ax.set_title('DFT(PBE) ecut conv.',fontsize=label_fs)
  ax.set_ylabel('total energy (Ry)',fontsize=label_fs)
  sel = (df['func'] == 'pbe')
  ax.plot(df.loc[sel,'ecut'],df.loc[sel,'E'])
  energy_ylim = ax.get_ylim()
  ax.axvline( df.loc[df['func']=='hse','ecut'].iloc[0] ,c='k',ls='--' )

  # magnetization vs. ecut
  mag_ylim = (4.,10.)
  iplot = 3
  ax = fig.add_subplot(2,2,iplot)
  ax.set_xlabel('ecut (Ry)',fontsize=label_fs)
  ax.set_ylabel('total magnetization (bohr mag.)',fontsize=label_fs)
  sel = (df['func'] == 'pbe')
  ax.plot(df.loc[sel,'ecut'],df.loc[sel,'abs_mag'])
  ax.set_ylim(mag_ylim)
  ax.axvline( df.loc[df['func']=='hse','ecut'].iloc[0] ,c='k',ls='--' )

  # magnetization vs. exx
  iplot = 4
  ax = fig.add_subplot(2,2,iplot)
  ax.set_xlabel('exact exchange fraction',fontsize=label_fs)
  ax.set_ylabel('total magnetization (bohr mag.)',fontsize=label_fs)
  sel = (df['func'] == 'hse')
  mydf = df[sel].sort_values('exx')
  ax.plot(mydf['exx'],mydf['abs_mag'])
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position("right")
  ax.set_ylim(mag_ylim)

  # total energy vs. exx
  iplot = 2
  ax = fig.add_subplot(2,2,iplot)
  ax.set_title('DFT(HSE) mag. vs. exx',fontsize=label_fs)
  ax.set_ylabel('total energy (Ry)',fontsize=label_fs)
  sel = (df['func'] == 'hse')
  mydf = df[sel].sort_values('exx')
  ax.plot(mydf['exx'],mydf['E'])
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position("right")
  ax.set_ylim(energy_ylim)

  plt.show()
# end __main__

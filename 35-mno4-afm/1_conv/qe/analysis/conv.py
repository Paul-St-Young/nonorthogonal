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

def plot1_mag_vs_exx(df):
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
# end def plot1_mag_vs_exx

def plot2_mag_vs_u(df):
  label_fs = 14
  fig = plt.figure()

  ax  = fig.add_subplot(2,2,1)
  ax.plot(df.loc[sel,'myu'],df.loc[sel,'E'])
  ax.set_ylabel('total energy (Ry)',fontsize=label_fs)

  ax  = fig.add_subplot(2,2,3)
  ax.set_xlabel('Hubbard U (eV)')
  ax.set_ylabel('absolute magnetization (a.u.)',fontsize=label_fs)
  ax.plot(df.loc[sel,'myu'],df.loc[sel,'abs_mag'])

  ax  = fig.add_subplot(2,2,2)
  ax.plot(df.loc[sel,'myu'],df.loc[sel,'pressure'])
  ax.set_ylabel('pressure (kbar)',fontsize=label_fs)
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position("right")

  ax  = fig.add_subplot(2,2,4)
  ax.plot(df.loc[sel,'myu'],df.loc[sel,'tot_mag'])
  ax.set_xlabel('Hubbard U (eV)')
  ax.set_ylabel('total magnetization (a.u.)',fontsize=label_fs)
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position("right")

  plt.show()
# end def plot2_mag_vs_u

def myu_from_input(inp_dict):
  myu  = float(inp_dict['system']['hubbard_u(1)'])
  return myu

if __name__ == '__main__':
  df0 = pd.read_json('../mno4_qe_ecut_exx.json')

  # get func, ecut and exx
  settings = df0['input'].apply(func_ecut_exx_from_input)
  df = pd.concat([settings,df0],axis=1)

  plot1_mag_vs_exx(df)

  # get hubbard U
  sel = (df['func'] == 'lda')
  df['myu'] = np.nan
  df.loc[sel,'myu'] = df.loc[sel,'input'].apply(myu_from_input)

  plot2_mag_vs_u(df)


# end __main__

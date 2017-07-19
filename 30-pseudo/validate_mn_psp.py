#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':

  ls_map = {'hgh':'-','jtk_s':':','jtk_p':'--'}
  cmap   = {'hgh':'k','jtk_s':'r','jtk_p':'b'}

  df0 = pd.read_json('3_mno_core_ref/core-hgh.json')
  df0['psp'] = 'hgh'
  df1 = pd.read_json('4_mno_jtk_slocal/jtk_s-bfd.json')
  df1['psp'] = 'jtk_s'
  df2 = pd.read_json('5_mno_jtk_plocal/jtk_p-bfd.json')
  df2['psp'] = 'jtk_p'

  df = pd.concat([df0,df1,df2]).reset_index(drop=True)
  fig,ax = plt.subplots(1,1)
  ax.set_xlabel('volume (a.u.)',fontsize=16)
  ax.set_ylabel('DFT(PBE) energy relative to volume 294 a.u. (ha)',fontsize=16)

  for psp in df.psp.unique():
    mydf = df[df.psp == psp].sort_values('volume')
    myx = mydf['volume'].values
    myy = mydf['E'].values

    # !!!! hack to get reference entry
    ref_sel = (mydf['infile_name'] == 'pbe-100-scf.in')
    myy0 = mydf.loc[ref_sel,'E'].values[0]

    ax.plot(myx,myy-myy0,label=psp,ls=ls_map[psp],c=cmap[psp],marker='.')
  # end for
  ax.legend(loc='upper left')
  fig.tight_layout()
  fig.savefig('mn_psp_validation.png',dpi=300)
  plt.show()

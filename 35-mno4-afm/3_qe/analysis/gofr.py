#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np

def gofr_from_file(fname,nequil,name):
  dmc_df = path_value_dataframe(fname,nequil)
  sel = dmc_df['h5path'].apply(lambda x:x.startswith(name))

  # calculate mean
  gr  = dmc_df.loc[sel,'value_mean'].values[0]
  gre = dmc_df.loc[sel,'value_error'].values[0]

  return gr,gre
# end def

def get_exx_gofr(rundir,nequil,name,exx_list=[0,25,50,75,100]):
  mean_list = []
  error_list = []
  for exx in exx_list:
    folder = os.path.join(rundir,'exx','hse-%d-ecut320'%exx,'dmc')
    flist  = sp.check_output(['find',folder,'-path','*.stat.h5']).split('\n')[:-1]
    fname = flist[-1]

    mymean,myerror = gofr_from_file(fname,nequil,name)

    mean_list.append(mymean)
    error_list.append(myerror)
  # end for
  return exx_list,mean_list,error_list
  #return np.array([exx_list,mean_list,error_list]).T
# end def

if __name__ == '__main__':
  import sys
  sys.path.insert(0,'../../../utils')
  from gather_stat import path_value_dataframe

  rundir = '../cont_scan'
  nequil = 200 # blocks
  name   = 'gofr_e_0_1'
  rmax   = 2.522759024

  for name in ['gofr_e_0_0','gofr_e_0_1','gofr_e_1_1']:
    exx_list,mean_list,error_list = get_exx_gofr(rundir,nequil,name)

    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel('r (bohr)',fontsize=16)
    ax.set_ylabel(name,fontsize=16)
    for iexx in range(len(exx_list)):
      myy = mean_list[iexx]
      myye=error_list[iexx]
      myx = np.linspace(0.0,rmax,len(myy))
      ax.errorbar(myx,myy,yerr=myye,label='exx=%d'%exx_list[iexx])
    # end for iexx
    ax.legend(loc=0)
    fig.tight_layout()
    fig.savefig('exx_%s.png'%name,dpi=320)
    plt.show()
  # end for

# end __main__

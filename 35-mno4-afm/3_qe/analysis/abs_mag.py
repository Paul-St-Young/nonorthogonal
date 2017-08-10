#!/usr/bin/env python
import os
import subprocess as sp
import numpy as np

def mabs_from_file(fname,nequil):
  dmc_df = path_value_dataframe(fname,nequil)
  sel_up = dmc_df['h5path'].apply(lambda x:x.startswith('spin_density') and ('/u/' in x))
  sel_dn = dmc_df['h5path'].apply(lambda x:x.startswith('spin_density') and ('/d/' in x))

  # calculate mean
  rho_up = dmc_df.loc[sel_up,'value_mean'].values[0]
  rho_dn = dmc_df.loc[sel_dn,'value_mean'].values[0]
  mabs   = abs(rho_up-rho_dn).sum()

  # calculate error
  rho_upe = dmc_df.loc[sel_up,'value_error'].values[0]
  rho_dne = dmc_df.loc[sel_dn,'value_error'].values[0]
  rho_umde= np.sqrt( rho_upe**2.+rho_dne**2. )
  mabs_err= np.sqrt( (rho_umde**2.).sum() )
  return mabs,mabs_err
# end def

def fqmc_from_flist(flist,iqmc):
  qstr = 's'+str(iqmc).zfill(3)
  for fname in flist:
    if qstr in fname:
      return fname
  raise RuntimeError('series %d not found in %s'%(iqmc,str(flist)) )
# end def

def get_exx_mabs(rundir,nequil,exx_list=[0,25,50,75,100],iqmc=2):
  mabs_list = []
  mabs_err_list = []
  for exx in exx_list:
    folder = os.path.join(rundir,'exx','hse-%d-ecut320'%exx,'dmc')
    flist  = sp.check_output(['find',folder,'-path','*.stat.h5']).split('\n')[:-1]
    fname = fqmc_from_flist(flist,iqmc)

    mabs,mabs_err = mabs_from_file(fname,nequil)

    mabs_list.append(mabs)
    mabs_err_list.append( mabs_err )
  # end for
  return np.array([exx_list,mabs_list,mabs_err_list]).T
# end def

def get_ldau_mabs(rundir,nequil,u_list=[0.5,2.0,4.0,5.0,6.0,8.0,10.0],iqmc=2):
  myu_list  = []
  mabs_list = []
  mabs_err_list = []
  for myu in u_list:
    folder = os.path.join(rundir,'ldau','u%3.2f'%myu,'dmc')
    flist  = sp.check_output(['find',folder,'-path','*.stat.h5']).split('\n')[:-1]
    fname = fqmc_from_flist(flist,iqmc)

    try:
      mabs,mabs_err = mabs_from_file(fname,nequil)
    except:
      continue

    myu_list.append(myu)
    mabs_list.append(mabs)
    mabs_err_list.append( mabs_err )
  # end for
  return np.array([myu_list,mabs_list,mabs_err_list]).T
# end def

def get_noci_mabs(rundir,nequil,ndet_list=[1,25,50,100,150],iqmc=2):
  my_ndet_list = []
  mabs_list = []
  mabs_err_list = []
  for idet in ndet_list:
    folder = os.path.join(rundir,'det%d'%idet)
    flist  = sp.check_output(['find',folder,'-path','*.stat.h5']).split('\n')[:-1]
    try:
      fname = fqmc_from_flist(flist,iqmc)
    except:
      continue

    mabs,mabs_err = mabs_from_file(fname,nequil)

    my_ndet_list.append(idet)
    mabs_list.append(mabs)
    mabs_err_list.append( mabs_err )
  # end for
  return np.array([my_ndet_list,mabs_list,mabs_err_list]).T
# def get_noci_mabs

if __name__ == '__main__':
  import sys
  sys.path.insert(0,'../../../utils')
  from gather_stat import path_value_dataframe

  fname = 'mabs_ldau.dat'
  if not os.path.isfile(fname):
    rundir = '../best_scan'
    nequil = 20 # blocks
    mabs_arr = get_ldau_mabs(rundir,nequil,iqmc=0)
    np.savetxt('vmc_'+fname,mabs_arr)
    nequil = 200 # blocks
    mabs_arr = get_ldau_mabs(rundir,nequil)
    np.savetxt(fname,mabs_arr)
  # end if

  exx_fname = 'mabs_exx.dat'
  if not os.path.isfile(exx_fname):
    rundir = '../best_scan'
    nequil = 200 # blocks
    hse_exx_mabs_arr = get_exx_mabs(rundir,nequil)
    np.savetxt(exx_fname,hse_exx_mabs_arr)
  # end if

  noci_fname = 'mabs_noci.dat'
  if not os.path.isfile(noci_fname):
    rundir = '../../../38-pyscf-basis-gs/2_dz_gs40/e-mdet/dmc'
    nequil = 20 # blocks
    mabs_arr = get_noci_mabs(rundir,nequil)
    np.savetxt(noci_fname,mabs_arr)
    mabs_arr = get_noci_mabs(rundir,nequil,iqmc=0)
    np.savetxt('vmc_'+noci_fname,mabs_arr)
  # end if

# end __main__
